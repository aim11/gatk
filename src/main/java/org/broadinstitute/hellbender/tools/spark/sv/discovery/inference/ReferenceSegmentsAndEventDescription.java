package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import scala.Tuple2;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants.*;

/**
 * This struct contains two key pieces of information that provides interpretation of the event:
 * <p>
 *     Ordered list of reference segments on the event primary chromosome that
 *     are bounded by segmenting locations extracted above from {@link CpxVariantInducingAssemblyContig.Jump}'s.
 *
 *     The order in which these segments are stored in the list
 *     is the same as how they are tiled on the reference, and
 *     two neighboring segments always share a boundary base,
 *     e.g. {(chr1, 1, 100), (chr1, 100, 150), (chr1, 150, 180), ...}
 *     (this shared base is because the segmenting locations are extracted from jumps,
 *      but didn't keep track of the order of the originating jumps.)
 * </p>
 *
 * <p>
 *     Description of how the reference segments are arranged on the sample, including orientation changes.
 *     Each description must be one of the three:
 *     <ul>
 *         <li>
 *             a signed integer, indicating which reference segment is placed,
 *             and if negative, meaning the segment's orientation is flipped on the sample.
 *             The absolute value is 1 + index into the above segments.
 *         </li>
 *         <li>
 *             a string that conforms to the format by {@link SimpleInterval#toString()}
 *             for describing a string of sequence that could map to that particular location,
 *             which is disjoint from the region returned by {@link CpxVariantInducingAssemblyContig.BasicInfo#getRefRegionBoundedByAlphaAndOmega()}.
 *         </li>
 *         <li>
 *             a string literal of the form {@link #UNMAPPED_INSERTION%d},
 *             used for indicating a part of the assembly contig is uncovered by any (high quality) mappings,
 *             with %d indicating how many bases unmapped
 *         </li>
 *     </ul>
 * </p>
 *
 * Two Cpx are equal if they are the same in their
 * <ul>
 *     <li>
 *         affected range
 *     </li>
 *     <li>
 *         segments
 *     </li>
 *     <li>
 *         alt arrangements
 *     </li>
 *     <li>
 *         alt haplotype sequence
 *     </li>
 * </ul>
 */
@DefaultSerializer(ReferenceSegmentsAndEventDescription.Serializer.class)
final class ReferenceSegmentsAndEventDescription {

    // used for indicating a part of the assembly contig is uncovered by any (high quality) mappings.
    public static final String UNMAPPED_INSERTION = "UINS";

    private final SimpleInterval affectedRefRegion;
    private final List<SimpleInterval> referenceSegments;
    private final List<String> eventDescriptions;
    private final byte[] altSeq;

    /**
     * This is a special case where a contig has all of the middle alignment(s) mapped to some disjoint places
     * (different chromosome, or same chromosome but not in the region returned by
     * {@link CpxVariantInducingAssemblyContig.BasicInfo#getRefRegionBoundedByAlphaAndOmega()}, aka valid region),
     * AND the first and last alignment share a single base on their ref span
     * (hence we have only one jump location on the event primary chromosome).
     * Hence the {@code segmentingLocations} is a one single-base entry list.
     */
    ReferenceSegmentsAndEventDescription(final AssemblyContigWithFineTunedAlignments preprocessedTig,
                                         final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                         final List<AlignmentInterval> contigAlignments,
                                         final SimpleInterval singleBase) {

        if ( singleBase.size() != 1)
            throw new CpxVariantInterpreter.UnhandledCaseSeen(
                    "run into unseen case where only one reference segmenting location is found but its size is not 1:\n"
                            + contigAlignments.toString());

        final SimpleInterval refRegionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();
        final boolean allMiddleAlignmentsDisjointFromAlphaOmega =
                contigAlignments
                        .subList(1, contigAlignments.size() - 1).stream()
                        .allMatch(ai -> CpxVariantInducingAssemblyContig.alignmentIsDisjointFromAlphaOmega(ai.referenceSpan, refRegionBoundedByAlphaAndOmega));
        if ( ! allMiddleAlignmentsDisjointFromAlphaOmega )
            throw new CpxVariantInterpreter.UnhandledCaseSeen("run into unseen case where only one reference segmenting location is found" +
                    " but some middle alignments are overlapping alpha-omega region:\t" + refRegionBoundedByAlphaAndOmega + "\n"
                    + contigAlignments.toString());

        affectedRefRegion = singleBase;
        referenceSegments = new ArrayList<>(Collections.singletonList(singleBase));

        eventDescriptions =
                contigAlignments
                        .subList(1, contigAlignments.size() - 1).stream()
                        .map(ai -> {
                            if (basicInfo.forwardStrandRep)
                                return ai.forwardStrand ? ai.referenceSpan.toString() : "-"+ai.referenceSpan.toString();
                            else
                                return ai.forwardStrand ? "-"+ai.referenceSpan.toString() : ai.referenceSpan.toString();
                        })
                        .collect(Collectors.toList());

        // the single base overlap from head and tail
        eventDescriptions.add(0, "1");
        eventDescriptions.add("1");

        altSeq = extractAltHaplotypeSeq(preprocessedTig, referenceSegments, basicInfo);
    }

    /**
     * For dealing with case where at least one middle alignments of the assembly contig are overlapping with
     * the region returned by {@link CpxVariantInducingAssemblyContig.BasicInfo#getRefRegionBoundedByAlphaAndOmega()}.
     * This is contrary to the case handled in
     * {@link #ReferenceSegmentsAndEventDescription(AssemblyContigWithFineTunedAlignments, CpxVariantInducingAssemblyContig.BasicInfo, List, SimpleInterval)}.
     */
    ReferenceSegmentsAndEventDescription(final CpxVariantInducingAssemblyContig cpxVariantInducingAssemblyContig) {

        final CpxVariantInducingAssemblyContig.BasicInfo basicInfo = cpxVariantInducingAssemblyContig.getBasicInfo();
        final List<AlignmentInterval> contigAlignments = cpxVariantInducingAssemblyContig.getPreprocessedTig().getSourceContig().alignmentIntervals;
        final List<CpxVariantInducingAssemblyContig.Jump> jumps = cpxVariantInducingAssemblyContig.getJumps();
        final List<SimpleInterval> segmentingLocations = cpxVariantInducingAssemblyContig.getEventPrimaryChromosomeSegmentingLocations();

        affectedRefRegion = getAffectedReferenceRegion(segmentingLocations);
        referenceSegments = extractRefSegments(basicInfo, segmentingLocations);
        eventDescriptions = makeInterpretation(basicInfo, contigAlignments, jumps, referenceSegments);

        altSeq = extractAltHaplotypeSeq(cpxVariantInducingAssemblyContig.getPreprocessedTig(), referenceSegments, basicInfo);
    }

    @VisibleForTesting
    static List<SimpleInterval> extractRefSegments(final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                                   final List<SimpleInterval> segmentingLocations) {

        final String eventPrimaryChromosome = basicInfo.eventPrimaryChromosome;

        final List<SimpleInterval> segments = new ArrayList<>(segmentingLocations.size() - 1);
        final Iterator<SimpleInterval> iterator = segmentingLocations.iterator();
        SimpleInterval leftBoundary = iterator.next();
        while (iterator.hasNext()) {
            final SimpleInterval rightBoundary = iterator.next();
            // there shouldn't be a segment constructed if two segmenting locations are adjacent to each other on the reference
            // this could happen when (in the simplest case), two alignments are separated by a mapped insertion (hence 3 total alignments),
            // and the two alignments' ref span are connected
            if (rightBoundary.getStart() - leftBoundary.getEnd() > 1) {
                segments.add(new SimpleInterval(eventPrimaryChromosome, leftBoundary.getStart(), rightBoundary.getStart()));
            }
            leftBoundary = rightBoundary;
        }
        return segments;
    }

    @VisibleForTesting
    static List<String> makeInterpretation(final CpxVariantInducingAssemblyContig.BasicInfo basicInfo,
                                           final List<AlignmentInterval> contigAlignments,
                                           final List<CpxVariantInducingAssemblyContig.Jump> jumps,
                                           final List<SimpleInterval> segments) {

        // using overlap with alignments ordered along the '+' strand representation of
        // the signaling contig to make sense of how the reference segments are ordered,
        // including orientations--using signs of the integers, on the sample;
        final List<AlignmentInterval> alignmentIntervalList =
                basicInfo.forwardStrandRep ? contigAlignments
                                           : ImmutableList.copyOf(contigAlignments).reverse();
        final Iterator<Integer> jumpGapSizeIterator =
                basicInfo.forwardStrandRep ? jumps.stream().map(jump ->jump.gapSize).iterator()
                                           : ImmutableList.copyOf(jumps).reverse().stream().map(jump ->jump.gapSize).iterator();
        Integer currentJumpGapSize = jumpGapSizeIterator.next();
        boolean jumpIsLast = false;

        final SimpleInterval regionBoundedByAlphaAndOmega = basicInfo.getRefRegionBoundedByAlphaAndOmega();

        final List<String> descriptions = new ArrayList<>(2*segments.size()); //ini. cap. a guess
        final List<Tuple2<SimpleInterval, Integer>> insertionsMappedToDisjointRegionsAndInsertionLocations = new ArrayList<>();
        for (final AlignmentInterval alignment : alignmentIntervalList) {

            if (CpxVariantInducingAssemblyContig.alignmentIsDisjointFromAlphaOmega(alignment.referenceSpan, regionBoundedByAlphaAndOmega)) {
                // disjoint alignment won't overlap any segments, so note down once where to insert, then move to next alignment
                final int indexABS = descriptions.size();
                insertionsMappedToDisjointRegionsAndInsertionLocations.add(new Tuple2<>(alignment.referenceSpan,
                        basicInfo.forwardStrandRep == alignment.forwardStrand ? indexABS : -1*indexABS));
                if ( currentJumpGapSize > 0 )
                    descriptions.add(UNMAPPED_INSERTION + "-" + currentJumpGapSize);
            } else {
                // depending on the representation and the current alignment's orientation, traverse segments in different order
                final int start, stop, step;
                if (basicInfo.forwardStrandRep == alignment.forwardStrand) {
                    start = 0;
                    stop = segments.size();
                    step = 1;
                } else {
                    start = segments.size()-1;
                    stop = -1;
                    step = -1;
                }
                // N*M overlaps
                for ( int i = start; i != stop; i += step ) {
                    final SimpleInterval currentSegment = segments.get(i);
                    // if current segment is contained in current alignment, note it down
                    if ( alignment.referenceSpan.contains(currentSegment) ) {
                        if (basicInfo.forwardStrandRep) // +1 below on i for 1-based description, no magic
                            descriptions.add( String.valueOf((alignment.forwardStrand ? 1 : -1) * (i+1)) );
                        else
                            descriptions.add( String.valueOf((alignment.forwardStrand ? -1 : 1) * (i+1)) );
                    }

                }
                // if the current alignment is associated with a gapped jump,
                // we need to signal that an unmapped insertion is present,
                // but only under two cases:
                //  1) no more segments to explore
                //  2) the next segment IS NOT contained in the current alignment's ref span
                if ( currentJumpGapSize > 0 && !jumpIsLast){
                    descriptions.add(UNMAPPED_INSERTION + "-" + currentJumpGapSize);
                }
            }

            if (jumpGapSizeIterator.hasNext()) // last alignment has no leaving jump so need to guard against that
                currentJumpGapSize = jumpGapSizeIterator.next();
            else
                jumpIsLast = true;
        }

        // post-hoc treatment of insertions that map to disjoint regions
        // go in reverse order because inserting into "descriptions", and insertion invalidates indices/iterators
        final ImmutableList<Tuple2<SimpleInterval, Integer>> reverse =
                ImmutableList.copyOf(insertionsMappedToDisjointRegionsAndInsertionLocations).reverse();
        for (final Tuple2<SimpleInterval, Integer> pair : reverse){
            final int index = pair._2;
            if (index > 0) {
                descriptions.add(pair._2, pair._1.toString());
            } else {
                descriptions.add(-1*pair._2, "-"+pair._1.toString());
            }
        }
        return descriptions;
    }

    /**
     * Extract alt haplotype sequence from the {@code contigWithFineTunedAlignments} to accompany the interpreted events.
     */
    @VisibleForTesting
    static byte[] extractAltHaplotypeSeq(final AssemblyContigWithFineTunedAlignments tigWithInsMappings,
                                         final List<SimpleInterval> segments,
                                         final CpxVariantInducingAssemblyContig.BasicInfo basicInfo) {

        final AlignmentInterval head = tigWithInsMappings.getSourceContig().getHeadAlignment();
        final AlignmentInterval tail = tigWithInsMappings.getSourceContig().getTailAlignment();
        if (head == null || tail == null)
            throw new GATKException("Head or tail alignment is null from contig:\n" + tigWithInsMappings.getSourceContig().toString());

        if (segments.isEmpty()) { // case where middle alignments all map to disjoint locations
            final int start = head.endInAssembledContig;
            final int end = tail.startInAssembledContig;
            final byte[] altSeq = Arrays.copyOfRange(tigWithInsMappings.getSourceContig().contigSequence, start - 1, end);
            if ( ! basicInfo.forwardStrandRep ) {
                SequenceUtil.reverseComplement(altSeq);
            }
            return altSeq;
        }

        final SimpleInterval firstSegment, lastSegment;
        if ( basicInfo.forwardStrandRep ) {
            firstSegment = segments.get(0);
            lastSegment = segments.get(segments.size() - 1);
        } else {
            firstSegment = segments.get(segments.size() - 1);
            lastSegment = segments.get(0);
        }

        final int start, end;

        if ( !firstSegment.overlaps(head.referenceSpan) ) {
            // if first segment doesn't overlap with head alignment,
            // it must be the case that the base (and possibly following bases) immediately after (or before if reverse strand) the head alignment's ref span is deleted
            final boolean firstSegmentNeighborsHeadAlignment = basicInfo.forwardStrandRep ? (firstSegment.getStart() - head.referenceSpan.getEnd() == 1)
                    : (head.referenceSpan.getStart() - firstSegment.getEnd() == 1);
            if ( ! firstSegmentNeighborsHeadAlignment )
                throw new CpxVariantInterpreter.UnhandledCaseSeen("1st segment is not overlapping with head alignment but it is not immediately before/after the head alignment either\n"
                        + tigWithInsMappings.toString());
            start = head.endInAssembledContig;
        } else {
            final SimpleInterval intersect = firstSegment.intersect(head.referenceSpan);
            if (intersect.size() == 1) {
                start = head.endInAssembledContig;
            } else {
                start = head.readIntervalAlignedToRefSpan(intersect)._1;
            }
        }

        if ( !lastSegment.overlaps(tail.referenceSpan) ) {
            final boolean expectedCase = basicInfo.forwardStrandRep ? (tail.referenceSpan.getStart() - lastSegment.getEnd() == 1)
                    : (lastSegment.getStart() - tail.referenceSpan.getEnd() == 1);
            if ( ! expectedCase )
                throw new CpxVariantInterpreter.UnhandledCaseSeen(tigWithInsMappings.toString());
            end = tail.startInAssembledContig;
        } else {
            final SimpleInterval intersect = lastSegment.intersect(tail.referenceSpan);
            if (intersect.size() == 1) {
                end = tail.startInAssembledContig;
            } else {
                end = tail.readIntervalAlignedToRefSpan(intersect)._2;
            }
        }

        // note from 1-based inclusive coordinate to C-style coordinate
        final byte[] altSeq = Arrays.copyOfRange(tigWithInsMappings.getSourceContig().contigSequence, start - 1, end);
        if ( ! basicInfo.forwardStrandRep ) {
            SequenceUtil.reverseComplement(altSeq);
        }

        return altSeq;
    }

    @VisibleForTesting
    static SimpleInterval getAffectedReferenceRegion(final List<SimpleInterval> eventPrimaryChromosomeSegmentingLocations) {
        final int start = eventPrimaryChromosomeSegmentingLocations.get(0).getStart();
        final int end = eventPrimaryChromosomeSegmentingLocations.get(eventPrimaryChromosomeSegmentingLocations.size() - 1).getEnd();
        return new SimpleInterval(eventPrimaryChromosomeSegmentingLocations.get(0).getContig(), start, end);
    }

    @VisibleForTesting
    VariantContextBuilder toVariantContext(final ReferenceMultiSource reference) throws IOException {

        final CpxVariantType cpxVariant = new CpxVariantType(affectedRefRegion, typeSpecificExtraAttributes());
        final SimpleInterval pos = new SimpleInterval(affectedRefRegion.getContig(), affectedRefRegion.getStart(), affectedRefRegion.getStart());

        final VariantContextBuilder vcBuilder = new VariantContextBuilder()
                .chr(affectedRefRegion.getContig()).start(affectedRefRegion.getStart()).stop(affectedRefRegion.getEnd())
                .alleles(AnnotatedVariantProducer.produceAlleles(pos, reference, cpxVariant))
                .id(cpxVariant.getInternalVariantId())
                .attribute(SVTYPE, cpxVariant.toString())
                .attribute(VCFConstants.END_KEY, affectedRefRegion.getEnd())
                .attribute(SVLEN, cpxVariant.getSVLength())
                .attribute(SEQ_ALT_HAPLOTYPE, new String(altSeq));

        cpxVariant.getTypeSpecificAttributes().forEach(vcBuilder::attribute);
        return vcBuilder;
    }

    private Map<String, String> typeSpecificExtraAttributes() {
        final Map<String, String> attributes = new HashMap<>();
        attributes.put(CPX_EVENT_ALT_ARRANGEMENTS, String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, eventDescriptions));
        if ( ! referenceSegments.isEmpty() ) {
            attributes.put(CPX_SV_REF_SEGMENTS,
                    String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR,
                            referenceSegments.stream().map(SimpleInterval::toString).collect(Collectors.toList())));
        }
        return attributes;
    }

    // =================================================================================================================

    List<SimpleInterval> getReferenceSegments() {
        return referenceSegments;
    }
    List<String> getEventDescriptions() {
        return eventDescriptions;
    }
    byte[] getAltSeq() {
        return altSeq;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder("ReferenceSegmentsAndEventDescription{");
        sb.append("affectedRefRegion=").append(affectedRefRegion);
        sb.append(", referenceSegments=").append(referenceSegments);
        sb.append(", eventDescriptions=").append(eventDescriptions);
        sb.append(", altSeq=").append(Arrays.toString(altSeq));
        sb.append('}');
        return sb.toString();
    }

    private ReferenceSegmentsAndEventDescription(final Kryo kryo, final Input input) {

        final String refRegionChr = input.readString();
        final int refRegionStart = input.readInt();
        final int refRegionEnd = input.readInt();
        affectedRefRegion = new SimpleInterval(refRegionChr, refRegionStart, refRegionEnd);

        final int numSegments = input.readInt();
        referenceSegments = new ArrayList<>(numSegments);
        for (int i = 0; i < numSegments; ++i)
            referenceSegments.add(kryo.readObject(input, SimpleInterval.class));

        final int numDescriptions = input.readInt();
        eventDescriptions = new ArrayList<>(numDescriptions);
        for (int i = 0; i < numDescriptions; ++i)
            eventDescriptions.add(input.readString());

        altSeq = new byte[input.readInt()];
        input.read(altSeq);
    }

    public void serialize(final Kryo kryo, final Output output) {
        output.writeString(affectedRefRegion.getContig());
        output.writeInt(affectedRefRegion.getStart());
        output.writeInt(affectedRefRegion.getEnd());

        output.writeInt(referenceSegments.size());
        for (final SimpleInterval segment : referenceSegments)
            kryo.writeObject(output, segment);

        output.writeInt(eventDescriptions.size());
        for (final String description: eventDescriptions)
            output.writeString(description);

        output.writeInt(altSeq.length);
        output.write(altSeq);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<ReferenceSegmentsAndEventDescription> {
        @Override
        public void write(final Kryo kryo, final Output output, final ReferenceSegmentsAndEventDescription x) {
            x.serialize(kryo, output);
        }

        @Override
        public ReferenceSegmentsAndEventDescription read(final Kryo kryo, final Input input,
                                                         final Class<ReferenceSegmentsAndEventDescription> clazz) {
            return new ReferenceSegmentsAndEventDescription(kryo, input);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final ReferenceSegmentsAndEventDescription that = (ReferenceSegmentsAndEventDescription) o;

        if (!referenceSegments.equals(that.referenceSegments)) return false;
        if (!affectedRefRegion.equals(that.affectedRefRegion)) return false;
        if (!eventDescriptions.equals(that.eventDescriptions)) return false;
        return Arrays.equals(altSeq, that.altSeq);
    }

    @Override
    public int hashCode() {
        int result = referenceSegments.hashCode();
        result = 31 * result + affectedRefRegion.hashCode();
        result = 31 * result + eventDescriptions.hashCode();
        result = 31 * result + Arrays.hashCode(altSeq);
        return result;
    }
}
