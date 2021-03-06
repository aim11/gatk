package org.broadinstitute.hellbender.tools.funcotator.dataSources;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.funcotator.Funcotation;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.LocatableXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.xsv.SimpleKeyXsvFuncotationFactory;
import org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer;
import org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable.XsvTableFeature;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A {@link Funcotation} to hold data from simple tabular data.
 * This class will be used any time an annotation is created by reading data from any data type
 * that can be expressed as a row in a table / database (such as via {@link SimpleKeyXsvFuncotationFactory},
 * {@link LocatableXsvFuncotationFactory}, and
 * {@link org.broadinstitute.hellbender.tools.funcotator.dataSources.cosmic.CosmicFuncotationFactory}).
 * Created by jonn on 11/28/17.
 */
public class TableFuncotation implements Funcotation {

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    /**
     * Names of the fields in this {@link TableFuncotation}.
     */
    private LinkedHashMap<String, String> fieldMap;

    //==================================================================================================================
    // Constructors:

    public TableFuncotation(final List<String> fieldNames, final List<String> fieldValues ) {
        if ( fieldNames.size() != fieldValues.size() ) {
            throw new UserException.BadInput("Field names and Field values are of different lengths!  This must not be!");
        }

        fieldMap = new LinkedHashMap<>(fieldNames.size());
        for ( int i = 0; i < fieldNames.size() ; ++i ) {
            fieldMap.put(fieldNames.get(i), fieldValues.get(i));
        }
    }

    public TableFuncotation(final XsvTableFeature xsvTableFeature) {

        final List<String> keys = xsvTableFeature.getHeaderWithoutLocationColuns();
        final List<String> values = xsvTableFeature.getValuesWithoutLocationColumns();

        fieldMap = new LinkedHashMap<>(keys.size());
        for ( int i = 0; i < keys.size() ; ++i ) {
            fieldMap.put(keys.get(i), values.get(i));
        }
    }

    public TableFuncotation(final LinkedHashMap<String, String> fieldMap ) {
        this.fieldMap = fieldMap;
    }

    public TableFuncotation(final TableFuncotation that ) {
        this.fieldMap = new LinkedHashMap<>( that.fieldMap );
    }

    //==================================================================================================================
    // Override Methods:

    @Override
    public void setFieldSerializationOverrideValue(final String fieldName, final String overrideValue) {
        if ( !fieldMap.containsKey(fieldName) ) {
            throw new GATKException("Attempted to override a field that is not contained in this TableFuncotation: "
                    + fieldName + " is not one of [" + String.join(",", fieldMap.keySet()) + "]");
        }
        fieldMap.put(fieldName, overrideValue);
    }

    @Override
    public String serializeToVcfString() {
        return fieldMap.values().stream()
                .map(f -> (f == null ? "" : f))
                .collect(Collectors.joining(VcfOutputRenderer.FIELD_DELIMITER));
    }

    @Override
    public LinkedHashSet<String> getFieldNames() {
        return new LinkedHashSet<>(fieldMap.keySet());
    }

    @Override
    public String getField(final String fieldName) {
        if ( fieldMap.containsKey(fieldName) ) {
            return fieldMap.get(fieldName);
        }
        else {
            throw new GATKException(this.getClass().getSimpleName() + ": Does not contain field: " + fieldName);
        }
    }

    @Override
    public boolean equals(final Object o) {
        if ( this == o ) return true;
        if ( o == null || getClass() != o.getClass() ) return false;

        final TableFuncotation that = (TableFuncotation) o;

        return fieldMap != null ? fieldMap.equals(that.fieldMap) : that.fieldMap == null;
    }

    @Override
    public int hashCode() {
        return fieldMap != null ? fieldMap.hashCode() : 0;
    }

    @Override
    public String toString() {
        return "TableFuncotation{" +
                "fieldMap={" + fieldMap.keySet().stream().map(k -> k + ":" + fieldMap.get(k)).collect(Collectors.joining(" , ")) + '}' +
                '}';
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    /**
     * Get the value in this {@link TableFuncotation} corresponding to the given key.
     * If the key is not in this {@link TableFuncotation}, returns {@code null}.
     * @param key Key to get from this {@link TableFuncotation}.
     * @return The value corresponding to the given key or {@code null}.
     */
    public String get(final String key) {
        return fieldMap.get(key);
    }

    /**
     * @return The {@link Set} of field names in this {@link TableFuncotation}.
     */
    public Set<String> keySet() {
        return fieldMap.keySet();
    }

    /**
     * @return The {@link Collection} of field values in this {@link TableFuncotation}.
     */
    public Collection<String> values() {
        return fieldMap.values();
    }

    /**
     * @return The number of key-value pairs in this {@link TableFuncotation}.
     */
    public int size() {
        return fieldMap.size();
    }

    //==================================================================================================================
    // Helper Data Types:
}
