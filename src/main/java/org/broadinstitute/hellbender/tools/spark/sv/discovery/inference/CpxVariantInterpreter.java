package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryInputData;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.ContigAlignmentsModifier;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.CONTIG_NAMES;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.CTG_GOOD_NONCANONICAL_MAPPING;
import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.TOTAL_MAPPINGS;

/**
 * This deals with the special case where a contig has multiple (> 2) alignments
 * and seemingly has the complete alt haplotype assembled.
 * See criteria in {@link AssemblyContigWithFineTunedAlignments#hasIncompletePictureFromMultipleAlignments()}.
 * For cases where the contig's alignment shows signature that the assembly doesn't
 * paint the full picture we could decide to emit all BND records,
 * but that could be dealt with later.
 */
public final class CpxVariantInterpreter {


    public static List<VariantContext> inferCpxVariant(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                                       final SvDiscoveryInputData svDiscoveryInputData) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryInputData.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryInputData.referenceSequenceDictionaryBroadcast;

        // almost every thing happens in this series of maps
        final JavaPairRDD<ReferenceSegmentsAndEventDescription, Iterable<CpxVariantInducingAssemblyContig>> interpretationAndAssemblyEvidence =
                assemblyContigs
                        .map(tig -> furtherPreprocess(tig, referenceSequenceDictionaryBroadcast.getValue()))
                        .map(tig -> new CpxVariantInducingAssemblyContig(tig, referenceSequenceDictionaryBroadcast.getValue()))
                        .mapToPair(tig -> new Tuple2<>(makeInterpretation(tig), tig))
                        .groupByKey(); // two contigs could give the same variant

        if (svDiscoveryInputData.discoverStageArgs.outputCpxResultsInHumanReadableFormat) {
            writeResultsForHumanConsumption(svDiscoveryInputData.outputPath, interpretationAndAssemblyEvidence);
        }

        return interpretationAndAssemblyEvidence.map(pair -> turnIntoVariantContext(pair, referenceBroadcast)).collect();
    }

    // =================================================================================================================

    /**
     * Essentially, this step is to de-overlap the alignments
     * (see {@link #deOverlapAlignments(List, SAMSequenceDictionary)})
     * because it would be very difficult to extract retracting jumps
     * (basically indicating homology, which is non-essential to event interpretation)
     * keeping track of them, and making sense of the event.
     *
     * @return the input contig with its alignments de-overlapped
     */
    private static AssemblyContigWithFineTunedAlignments furtherPreprocess(final AssemblyContigWithFineTunedAlignments contigWithFineTunedAlignments,
                                                                           final SAMSequenceDictionary refSequenceDictionary) {
        final AlignedContig sourceTig = contigWithFineTunedAlignments.getSourceContig();

        final List<AlignmentInterval> deOverlappedAlignmentConfiguration =
                deOverlapAlignments(sourceTig.alignmentIntervals, refSequenceDictionary);

        final AlignedContig contig = new AlignedContig(sourceTig.contigName, sourceTig.contigSequence,
                deOverlappedAlignmentConfiguration, sourceTig.hasEquallyGoodAlnConfigurations);

        return new AssemblyContigWithFineTunedAlignments(contig,
                contigWithFineTunedAlignments.getInsertionMappings(),
                contigWithFineTunedAlignments.getSAtagForGoodMappingToNonCanonicalChromosome());
    }

    @VisibleForTesting
    static List<AlignmentInterval> deOverlapAlignments(final List<AlignmentInterval> originalAlignments,
                                                       final SAMSequenceDictionary refSequenceDictionary) {
        final List<AlignmentInterval> result = new ArrayList<>(originalAlignments.size());
        final Iterator<AlignmentInterval> iterator = originalAlignments.iterator();
        AlignmentInterval one = iterator.next();
        while (iterator.hasNext()) {
            final AlignmentInterval two = iterator.next();
            // TODO: 11/5/17 an edge case is possible where the best configuration contains two alignments,
            //       one of which contains a large gap, and since the gap split happens after the configuration scoring,
            //       (that gap split happens after scoring is due to how MQ and AS are used in the scoring step, gap-split alignment cannot use originating alignment's values, but it takes time to recompute)
            //       one of the alignment from the gap split may be contained in the other original alignment, leading to problems;
            //       here we skip the alignment that is BEFORE the child alignment from the gap-split,
            //       IFF that alignment contains the child alignment in terms of their spans on the read/contig
            //       if you are concerned about the first child alignment from the same gapped alignment being skipped,
            //       don't worry, that is impossible because child alignments of the same gapped alignment cannot overlap on the read.
            if (two.alnModType.equals(ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT)) {
                final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);
                if (overlapOnRead >= two.getSizeOnRead())
                    continue;
            }
            final List<AlignmentInterval> deoverlapped = ContigAlignmentsModifier.removeOverlap(one, two, refSequenceDictionary);
            result.add(deoverlapped.get(0));
            one = deoverlapped.get(1);
        }
        result.add(one);
        return result;
    }

    @VisibleForTesting
    static ReferenceSegmentsAndEventDescription makeInterpretation(final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig) {

        if (cpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations().size() == 1) {
            return new ReferenceSegmentsAndEventDescription(
                    cpxVariantInducingAssemblyContig.getPreprocessedTig(),
                    cpxVariantInducingAssemblyContig.getBasicInfo(),
                    cpxVariantInducingAssemblyContig.getPreprocessedTig().getSourceContig().alignmentIntervals,
                    cpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations().get(0));
        } else {
            return new ReferenceSegmentsAndEventDescription(cpxVariantInducingAssemblyContig);
        }
    }

    @VisibleForTesting
    static VariantContext turnIntoVariantContext(final Tuple2<ReferenceSegmentsAndEventDescription, Iterable<CpxVariantInducingAssemblyContig>> pair,
                                                 final Broadcast<ReferenceMultiSource> referenceBroadcast)
            throws IOException {

        final VariantContextBuilder rawVariantContextBuilder = pair._1.toVariantContext(referenceBroadcast.getValue());

        final List<String> contigNames = new ArrayList<>();
        final List<String> mayContainNoInfo = new ArrayList<>(); // for storing AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME

        Utils.stream(pair._2)
                .sorted(Comparator.comparing(tig -> tig.getPreprocessedTig().getSourceContig().contigName))
                .forEach(evidenceTig -> {

                    contigNames.add(evidenceTig.getPreprocessedTig().getSourceContig().contigName);

                    final String saTagForGoodMappingToNonCanonicalChromosome =
                            evidenceTig.getPreprocessedTig().getSAtagForGoodMappingToNonCanonicalChromosome();
                    if ( ! saTagForGoodMappingToNonCanonicalChromosome
                            .equals(AssemblyContigWithFineTunedAlignments.NO_GOOD_MAPPING_TO_NON_CANONICAL_CHROMOSOME)) {
                        mayContainNoInfo.add(saTagForGoodMappingToNonCanonicalChromosome);
                    } else {
                        mayContainNoInfo.add(".");
                    }
                });

        final Map<String, String> attributeMap = new HashMap<>();

        attributeMap.put(TOTAL_MAPPINGS, String.valueOf(contigNames.size()));

        attributeMap.put(CONTIG_NAMES, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, contigNames));

        if ( ! mayContainNoInfo.equals( Collections.nCopies(contigNames.size(), ".") ) ) {
            attributeMap.put(CTG_GOOD_NONCANONICAL_MAPPING, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, mayContainNoInfo));
        }

        // TODO: 2/1/18 question: show we introduce new annotations for them?
//        attributeMap.put(GATKSVVCFConstants.HQ_MAPPINGS,       tigWithInsMappings.getSourceContig().alignmentIntervals.stream().mapToInt(ai -> ai.mapQual).anyMatch( mq -> mq < CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD ) ? "0" : "1");
//        attributeMap.put(GATKSVVCFConstants.MAPPING_QUALITIES, tigWithInsMappings.getSourceContig().alignmentIntervals.stream().map(ai -> String.valueOf(ai.mapQual)).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
//        attributeMap.put(GATKSVVCFConstants.ALIGN_LENGTHS,     tigWithInsMappings.getSourceContig().alignmentIntervals.stream().map(ai -> String.valueOf(ai.getSizeOnRead())).collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)));
//        attributeMap.put(GATKSVVCFConstants.MAX_ALIGN_LENGTH,  String.valueOf(tigWithInsMappings.getSourceContig().alignmentIntervals.stream().mapToInt(AlignmentInterval::getSizeOnRead).max().orElse(0)));

        // TODO: 12/11/17 integrate these with those that survived the alignment filtering step?
//        // known insertion mappings from filtered out alignments
//        if ( !tigWithInsMappings.getInsertionMappings().isEmpty() ) {
//            attributeMap.put(INSERTED_SEQUENCE_MAPPINGS,
//                    String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, tigWithInsMappings.getInsertionMappings()));
//        }

        attributeMap.forEach(rawVariantContextBuilder::attribute);
        return rawVariantContextBuilder.make();
    }

    // =================================================================================================================

    private static void writeResultsForHumanConsumption(final String outputPath,
                                                        final JavaPairRDD<ReferenceSegmentsAndEventDescription, Iterable<CpxVariantInducingAssemblyContig>> interpretationAndAssemblyEvidence) {
        try {
            // for easier view when debugging, will be taken out in the final commit.
            Files.write(Paths.get(Paths.get(outputPath).getParent().toAbsolutePath().toString() + "/cpxEvents.txt"),
                    () -> interpretationAndAssemblyEvidence
                            .flatMap(pair -> Utils.stream(pair._2).map( tig -> new Tuple2<>(tig, pair._1)).iterator())
                            .sortBy(pair  -> pair._1.getPreprocessedTig().getSourceContig().contigName, true, 1)
                            .map(pair -> {
                                final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig = pair._1;
                                final ReferenceSegmentsAndEventDescription referenceSegmentsAndEventDescription = pair._2;
                                String s = cpxVariantInducingAssemblyContig.toString() + "\n";
                                s += referenceSegmentsAndEventDescription.toString();
                                s += "\n";
                                return (CharSequence) s;
                            })
                            .collect().iterator());
        } catch (final IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile("Could not save filtering results to file", ioe);
        }
    }

    static final class UnhandledCaseSeen extends GATKException.ShouldNeverReachHereException {
        private static final long serialVersionUID = 0L;
        UnhandledCaseSeen( final String s ) {
            super(s);
        }

        UnhandledCaseSeen( final String s, final Throwable throwable ) {
            super(s, throwable);
        }

        UnhandledCaseSeen( final Throwable throwable) {this("Seeing unhandled case", throwable);}
    }
}
