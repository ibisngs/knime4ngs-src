package de.helmholtz_muenchen.ibis.utils.datatypes.file;


import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.codec.binary.Hex;
import org.apache.commons.codec.digest.DigestUtils;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;
import org.knime.core.data.DataValue;
import org.knime.core.data.StringValue;

@SuppressWarnings("serial")
public class FileCell extends DataCell implements FileValue, StringValue{
	   /** Type for File cells. */
    public static final DataType TYPE = DataType.getType(FileCell.class);

    /**
     * Returns the preferred value class of this cell implementation. This
     * method is called per reflection to determine which is the preferred
     * renderer, comparator, etc.
     *
     * @return FileValue.class;
     */
    public static final Class<? extends DataValue> getPreferredValueClass() {
        return FileValue.class;
    }

    private static final FileSerializer SERIALIZER = new FileSerializer();

    /**
     * Returns the factory to read/write DataCells of this class from/to a
     * DataInput/DataOutput. This method is called via reflection.
     *
     * @return a serializer for reading/writing cells of this kind
     * @see DataCell
     */
    public static final FileSerializer getCellSerializer() {
        return SERIALIZER;
    }

    private final String m_FileString;
    
    private String m_md5 = "";		//Has to be calculated
    
    private String m_fileExtension;
    private String m_Description;
    private String m_FileType;

    /**
     * Creates a new File Cell based on the given String value. This constructor
     * is used from the {@link FileCellFactory#create(String) create} method of
     * the accompanying {@link FileCellFactory}.
     *
     * @param str the String value to store
     * @throws NullPointerException if the given String value is
     *             <code>null</code>
     */
    FileCell(final String str) {
        if (str == null) {
            throw new NullPointerException("File value must not be null.");
        }
        m_FileString = str;
    }

    /** {@inheritDoc} */
    @Override
    public String getStringValue() {
        return m_FileString;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return getStringValue();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean equalsDataCell(final DataCell dc) {
        return m_FileString.equals(((FileCell)dc).m_FileString);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return m_FileString.hashCode();
    }

    /** Factory for (de-)serializing a FileCell. */
    private static class FileSerializer implements DataCellSerializer<FileCell> {
        /**
         * {@inheritDoc}
         */
        @Override
        public void serialize(final FileCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }

        /**
         * {@inheritDoc}
         */
        @Override
        public FileCell deserialize(final DataCellDataInput input)
                throws IOException {
            String s = input.readUTF();
            return new FileCell(s);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String getFilePath() {
        return m_FileString;
    }
    

    /**
     * Calculates the md5 checksum of FileCell
     * @param Overwrite Force to overwrite the old md5 by a new one
     */
    public void calculateMD5Value(boolean Overwrite){
    	if(m_md5.equals("")||Overwrite){	//Calculate only if md5_value has not been set yet or if Overwrite=true
    		m_md5 = computeMD5();
    	}else{
    		System.out.println("MD5 has already been calculated: "+m_md5);
    	}
    }

    /**
     * The md5 checksum
     * @return
     */
    public String getMD5Value(){
    	return m_md5;
    }
    
    /**
     * Private function to compute the md5 value
     * @return The computed MD5 value
     */
    private String computeMD5(){
    	
    	String md5value = "";
    	try{
    		//Read the file
	    	FileInputStream fis = new FileInputStream(new File(m_FileString));
	    	//Get the md5
	    	byte data[] = DigestUtils.md5(fis);
			char md5Chars[] = Hex.encodeHex(data);
			md5value = String.valueOf(md5Chars);

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	return md5value;
    }
    
    
    /**
     * Compares the original MD5 value vs. a currently computed one
     * @return True if MD5 value has no changed since last MD5 check
     */
    public boolean checkMD5Value(){
    	String tmpMD5 = computeMD5();
    	if(m_md5.equals(tmpMD5)){
    		if(m_md5.equals("")){
    			System.out.println("MD5 not set yet");
    		}else{
    			System.out.println("MD5 did not change.");
    		}
    		return true;
    	}else{
    		System.out.println("MD5 changed since last check !");
    		System.out.println("Old: "+m_md5);
    		System.out.println("New: "+tmpMD5);
    		return false;
    	}
    }
    
    /**
     * Sets the Meta Info for the FileCell 
     * @param ext	The file extension (.txt,.sam,...)
     * @param type	The file type (e.g. Samtools mpileup)
     * @param description A short file description
     */
    public void setMetaInfo(String ext,String type, String description){
    	m_fileExtension = ext;
    	m_FileType = type;
    	m_Description = description;
    }
    
    /**
     * Returns the file extension
     * @return File Extension
     */
    public String getFileExtension(){
    	return m_fileExtension;
    }
    
    /**
     * Returns the file description
     * @return File Description
     */
    public String getFileDescription(){
    	return m_Description;
    }
    
    /**
     * Returns the file type
     * @return
     */
    public String getFileType(){
    	return m_FileType;
    }
    
    
}
