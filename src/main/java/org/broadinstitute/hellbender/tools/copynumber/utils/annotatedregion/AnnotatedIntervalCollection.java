package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvLocatableTableCodec.getAndValidateConfigFileContents;

/**
 * Represents a collection of annotated regions.  The annotations do not need to be known ahead of time, if reading from a file.
 *
 *  This class supports reading xsv (tsv, csv) files with comments ("#") and SAM headers ("@").  The default is tsv.
 */
public class AnnotatedIntervalCollection {

    private static final Logger logger = LogManager.getLogger(AnnotatedIntervalCollection.class);

    public static final String ANNOTATED_REGION_DEFAULT_CONFIG_RESOURCE = "org/broadinstitute/hellbender/tools/copynumber/utils/annotatedregion/annotated_region_default.config";
    private final SAMFileHeader samFileHeader;

    /** Does not include the locatable fields. */
    private final List<String> annotations;

    /** Can be anything here.  These are comments for the AnnotatedIntervals. */
    private final List<String> comments;

    private final List<AnnotatedInterval> records;

    private AnnotatedIntervalCollection(final SAMFileHeader samFileHeader, final List<String> annotations, final List<String> comments,
                                        final List<AnnotatedInterval> records) {
        this.comments = comments;
        this.samFileHeader = samFileHeader;
        this.annotations = annotations;
        this.records = records;
    }

    /**
     *  Same as {@link #create(Path, Path, Set)} , but uses the default annotation
     *   region config file in the GATK.
     * @param input See {@link #create(Path, Path, Set)}
     * @param headersOfInterest See {@link #create(Path, Path, Set)}
     * @return See {@link #create(Path, Path, Set)}
     */
    public static AnnotatedIntervalCollection create(final Path input, final Set<String> headersOfInterest) {

        final String resourcePath = ANNOTATED_REGION_DEFAULT_CONFIG_RESOURCE;
        try {
            final File tmpResourceFile = Resource.getResourceContentsAsFile(resourcePath);
            return create(input, tmpResourceFile.toPath(), headersOfInterest);
        } catch (final IOException ioe) {
            throw new GATKException.ShouldNeverReachHereException("Could not read config file: " + resourcePath,
                    ioe);
        }
    }

    /**
     *  Returns a copy of the given annotated interval, but only with the annotations to preserve.
     *
     *  No checking is done to ensure that this will leave the copy with any annotations.
     *  No checking is done to ensure that any of the annotationsToPreserve are actually in the annotatedInterval.
     *
     * @param annotatedInterval
     * @param annotationsToPreserve
     * @return
     */
    private static AnnotatedInterval copyAnnotatedInterval(final AnnotatedInterval annotatedInterval, final List<String> annotationsToPreserve) {
        final SortedMap<String, String> copiedAnnotations = annotatedInterval.getAnnotations().entrySet().stream()
                .filter(e -> annotationsToPreserve.contains(e.getKey()))
                .collect(TreeMap::new, (map, e) -> map.put(e.getKey(), e.getValue()), (map, map2) -> map.putAll(map2));
        return new AnnotatedInterval(annotatedInterval.getInterval(), copiedAnnotations);
    }

    /**
     * Create a collection from components.
     *
     * @param regions regions to use in the resulting collection.  Never {@code null}.
     * @param samFileHeader SAMFileHeader to include in the collection.  Represents the sample(s)/references that were used for these segments.
     *                      {@code null} is allowed.
     * @param annotations List of annotations to preserve in the regions.  Never {@code null}.  These are the only annotations that will be written.
     * @param comments List of comments that need to be associated with this collection.  Use {@link Collections#emptyList()} when no comments needed.  Never {@code null}.
     * @return collection based on the inputs.  Never {@code null}.
     */
    public static AnnotatedIntervalCollection create(final List<AnnotatedInterval> regions,
                                                     final SAMFileHeader samFileHeader,
                                                     final List<String> annotations,
                                                     final List<String> comments) {

        Utils.nonNull(regions);
        Utils.nonNull(annotations);
        Utils.nonNull(comments);

        final List<AnnotatedInterval> updatedAnnotatedIntervals = regions.stream().map(r -> copyAnnotatedInterval(r, annotations))
                .collect(Collectors.toList());

        return new AnnotatedIntervalCollection(samFileHeader, annotations, comments, updatedAnnotatedIntervals);
    }

    //TODO: Update docs since this is not the primary way of creating a collection.  Move the detailed comments to the public method.
    /** Create a collection based on the contents of an input file and a given config file.  The config file must be the same as
     * is ingested by {@link XsvLocatableTableCodec}.
     *
     * @param input readable path to use for the xsv file.  Must be readable.  Never {@code null}.
     * @param inputConfigFile config file for specifying the format of the xsv file.  Must be readable.  Never {@code null}.
     * @param headersOfInterest Only preserve these headers.  These must be present in the input file.  This parameter should not include the locatable columns
     *                          defined by the config file, which are always preserved.
     *                          Use {@code null} to indicate "all headers are of interest".
     * @return never {@code null}
     */
    public static AnnotatedIntervalCollection create(final Path input, final Path inputConfigFile, final Set<String> headersOfInterest) {

        IOUtils.assertFileIsReadable(input);
        IOUtils.assertFileIsReadable(inputConfigFile);

        final XsvLocatableTableCodec codec = new XsvLocatableTableCodec(inputConfigFile);
        final List<AnnotatedInterval> regions = new ArrayList<>();

        if (codec.canDecode(input.toString())) {
            try (final InputStream fileInputStream = Files.newInputStream(input)) {

                // Lots of scaffolding to do reading here:
                final AsciiLineReaderIterator lineReaderIterator = new AsciiLineReaderIterator(AsciiLineReader.from(fileInputStream));
                final List<String> header = codec.readActualHeader(lineReaderIterator);
                warnAllHeadersOfInterestNotPresent(headersOfInterest, header);

                final List<String> annotationCols = codec.getHeaderWithoutLocationColumns();

                while (lineReaderIterator.hasNext()) {
                    final XsvTableFeature feature = codec.decode(lineReaderIterator.next());
                    if (feature == null) {
                        continue;
                    }

                    final List<String> featureValues = feature.getValuesWithoutLocationColumns();

                    final SortedMap<String, String> annotations = new TreeMap<>();
                    IntStream.range(0, annotationCols.size()).boxed()
                            .filter(i -> (headersOfInterest == null) || headersOfInterest.contains(annotationCols.get(i)))
                            .forEach(i -> annotations.put(annotationCols.get(i), featureValues.get(i)));

                    regions.add(new AnnotatedInterval(
                            new SimpleInterval(feature.getContig(), feature.getStart(), feature.getEnd()),
                            annotations));
                }

                return new AnnotatedIntervalCollection(codec.createSamFileHeader(), codec.getHeaderWithoutLocationColumns(),
                        codec.getComments(), regions);

            }
            catch ( final FileNotFoundException ex ) {
                throw new GATKException("Error - could not find file: " + input, ex);
            }
            catch ( final IOException ex ) {
                throw new GATKException("Error - IO problem with file " + input, ex);
            }
        }
        else {
            throw new UserException.BadInput("Could not parse xsv file.");
        }
    }

    private static void warnAllHeadersOfInterestNotPresent(final Set<String> headersOfInterest, final List<String> header) {
        if ((headersOfInterest != null) && !header.containsAll(headersOfInterest)) {
            final Set<String> unusedColumnsOfInterest = Sets.difference(new HashSet<>(headersOfInterest), new HashSet<>(header));
            if (unusedColumnsOfInterest.size() > 0) {
                final List<String> missingColumns = new ArrayList<>(unusedColumnsOfInterest);
                logger.warn("Some headers of interest specified by the user were not seen in input: " + StringUtils.join(missingColumns, ", "));
            }
        }
    }

    /**
     *  Write this collection to a file
     *
     *  Dev note:  This method will force the default xsv config file on the output.
     *
     * @param outputFile destination file, must be writable.
     */
    public void write(final File outputFile) {
        try (final AnnotatedIntervalWriter writer = new SimpleAnnotatedIntervalWriter(outputFile);){

            writer.writeHeader(createHeaderForWriter(annotations, samFileHeader, comments));
            getRecords().forEach(writer::add);
        } catch (final IOException ioe) {
            throw new GATKException("Error - could not write output file: " + outputFile.getAbsolutePath(), ioe);
        }
    }

    private static AnnotatedIntervalHeader createHeaderForWriter(final List<String> annotations, final SAMFileHeader samFileHeader, final List<String> comments) throws IOException {
        final File resourceFile = Resource.getResourceContentsAsFile(AnnotatedIntervalCollection.ANNOTATED_REGION_DEFAULT_CONFIG_RESOURCE);
        final Properties headerNameProperties = getAndValidateConfigFileContents(resourceFile.toPath());
        final String contigColumnName = headerNameProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_CONTIG_COLUMN_KEY);
        final String startColumnName = headerNameProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_START_COLUMN_KEY);
        final String endColumnName = headerNameProperties.getProperty(XsvLocatableTableCodec.CONFIG_FILE_END_COLUMN_KEY);

        XsvLocatableTableCodec.validateLocatableColumnName(contigColumnName);
        XsvLocatableTableCodec.validateLocatableColumnName(startColumnName);
        XsvLocatableTableCodec.validateLocatableColumnName(endColumnName);

        return new AnnotatedIntervalHeader(contigColumnName, startColumnName, endColumnName, annotations, samFileHeader, comments);
    }

    /** Can return {@code null} */
    public SAMFileHeader getSamFileHeader() {
        return samFileHeader;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    public List<String> getComments() {
        return comments;
    }

    public List<AnnotatedInterval> getRecords() {
        return records;
    }

    public int size() {
        return getRecords().size();
    }
}
