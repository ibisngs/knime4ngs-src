package de.helmholtz_muenchen.ibis.utils;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.InvalidSettingsException;

public class CompatibilityChecker {

	private String WarningMessage;
	private boolean Warning;
	
	
	public CompatibilityChecker(){
		
		WarningMessage = "NA";
		Warning 		= false;
		
	}
	
	/**
	 * Based on the inSpecs, the method checks for paired-end or single-end mapping 
	 * @param inSpecs
	 * @param port
	 * @return "paired-end" or "single-end"
	 * @throws InvalidSettingsException
	 */
	public String getReadType(DataTableSpec[] inSpecs, int port) throws InvalidSettingsException{
		
		int foundFastQCells = 0;
		
		for(int i=0;i<inSpecs[port].getNumColumns();i++){
			if(inSpecs[port].getColumnSpec(i).getType().toString().equals("FastQCell")){
				foundFastQCells++;
			}
		}
			
		if(foundFastQCells==1){
			return "single-end";
		}else if(foundFastQCells==2){
			return "paired-end";
		}else if (foundFastQCells==0){
			throw new InvalidSettingsException("No FastQCells in input table! Incompatible predecessor node?");
		}else{
			setWarningMessages("More than two FastQCells found. Using first two in paired-end mode.");
			return "paired-end";
		}
	}
	
	/**
	 * Returns a string of all warning messages
	 * @return
	 */
	public String getWarningMessages(){
		return WarningMessage;
	}
	
	/**
	 * Returns true if a warning message has been set. Use getWarningMessages() to retrieve it.
	 * @return
	 */
	public boolean getWarningStatus(){
		return Warning;
	}
	
	/**
	 * Appends the warning message to the global warning string
	 * @param warning
	 */
	private void setWarningMessages(String warning){
		if(getWarningMessages().equals("NA")){
			WarningMessage	= warning;
			Warning 		= true;
		}else{
			WarningMessage += "\n"+warning;
		}
	}
	
	public static boolean inputFileNotOk(String path) {
		return inputFileNotOk(path, true);
	}
	
	public static boolean inputFileNotOk(String path, boolean checkEmpty) {
		
		if(path.equals("") || path == null) return true;
		
		if(checkEmpty && Files.exists(Paths.get(path))) {
			boolean isEmpty = false;
			try {
				isEmpty = (Files.newBufferedReader(Paths.get(path),StandardCharsets.ISO_8859_1).readLine() == null);
			} catch (IOException e) {
				e.printStackTrace();
			}
			return isEmpty;
		}
		return Files.notExists(Paths.get(path));
	}
	
	public static boolean checkInputCellType(DataTableSpec inSpecs, String CellType) {
		return (getIndexCellType(inSpecs,CellType)>-1);
	}
	
	/**
	 * @deprecated
	 * @param inSpecs
	 * @param CellType
	 * @return
	 */
	public static int getIndexCellType(DataTableSpec inSpecs, String CellType) {
		int index = -1;
		
		for(int i = 0; i < inSpecs.getNumColumns(); i++) {
    		if(inSpecs.getColumnSpec(i).getType().toString().equals(CellType)) {
    			index = i;
    		}
    	}
		return index;
	}
	
	public static int getFirstIndexCellType(DataTableSpec inSpecs, String CellType) {
		int index = -1;
		
		for(int i = 0; i < inSpecs.getNumColumns(); i++) {
    		if(inSpecs.getColumnSpec(i).getType().toString().equals(CellType)) {
    			index = i;
    			break;
    		}
    	}
		return index;
	}
	
}
