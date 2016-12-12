
/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.helmholtz_muenchen.ibis.utils;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.InvalidSettingsException;

import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

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
		return inputFileNotOk(path, true, null);
	}
	
	public static boolean inputFileNotOk(String path, DataType type) {
		return inputFileNotOk(path, true, type);
	}
	
	public static boolean inputFileNotOk(String path, boolean checkEmpty) {
		return inputFileNotOk(path, checkEmpty, null);
	}
	
	public static boolean inputFileNotOk(String path, boolean checkEmpty, DataType type) {
		
		if( path == null || path.equals("")) return true;
		
		if(checkEmpty && Files.exists(Paths.get(path))) {
			boolean isEmpty = false;
			try {
				isEmpty = (Files.newBufferedReader(Paths.get(path),StandardCharsets.ISO_8859_1).readLine() == null);
			} catch (IOException e) {
				e.printStackTrace();
			}
			return isEmpty;
		}
		
		if(Files.notExists(Paths.get(path))) return true;
		
		if(type == null) {
			return false;
		} else {
			String my_type = type.toString();
			if(my_type.equals("FastACell")) {
				return !FileValidator.checkFastaFormat(path);
			}
			if(my_type.equals("FastQCell")) {
				return !FileValidator.checkFastqFormat(path);
			}
		}
		return false;
	}
	
	public static boolean checkInputCellType(DataTableSpec inSpecs, String CellType) {
		return (getFirstIndexCellType(inSpecs,CellType)>-1);
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
	
	/**
	 * Method checks for empty columns in data table and throws exception when no columns exist
	 * @param inSpecs
	 * @throws InvalidSettingsException
	 */
	
	public static void InputFileNoColumns(DataTableSpec[] inSpecs) throws InvalidSettingsException {
    	
    	for(int i = 0; i < inSpecs.length; i++ ){
    		
    		if(inSpecs[i].getNumColumns()==0){
    			throw new InvalidSettingsException ("One or more input tables have no columns");
    		}
    	}
	}
	
	/**
	 * Method checks for empty rows in data table and throws exception when no rows exists
	 * @param inData
	 * @throws InvalidSettingsException
	 */
	public static void InputFileNoRows(BufferedDataTable[] inData) throws InvalidSettingsException {
		
		for(int i = 0; i < inData.length; i++) {
			java.util.Iterator <DataRow> it = inData[i].iterator();
			
			if(it.hasNext()==false) {
					throw new InvalidSettingsException ("One or more input tables have no rows"); 
			}
		}
	
	}
	
	/**
	 * Method checks for missing values in input table and throws exception if empty string values exist
	 * @param inData
	 * @throws InvalidSettingsException
	 */
	private static void InputContainsEmptyCells(BufferedDataTable[] inData) throws InvalidSettingsException{
		for(int i = 0; i < inData.length; i++) {
			java.util.Iterator <DataRow> it = inData[i].iterator();
			
	        while(it.hasNext()){
	        	
	        	java.util.Iterator <DataCell> cellIT = it.next().iterator();
	        	
	        	while(cellIT.hasNext()){
	        		
	        		DataCell c = cellIT.next();
	        		if(c.toString().equals("") || c.isMissing()){
	        			throw new InvalidSettingsException ("One or more rows contain missing values!"); 
	        		}
	        	}
	        } 
			
			
			
		}
	}
	
	
	
	/**
	 * Method checks if inData is ok and ready for execution
	 * @param inData
	 * @throws InvalidSettingsException
	 */
	public static void inDataCheck(BufferedDataTable[] inData) throws InvalidSettingsException{
		//Check for empty table
		InputFileNoRows(inData);
		//Check missing values
		InputContainsEmptyCells(inData);
	}
	
	
}
