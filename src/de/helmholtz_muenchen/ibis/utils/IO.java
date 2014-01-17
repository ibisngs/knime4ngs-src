package de.helmholtz_muenchen.ibis.utils;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.URL;

import org.knime.base.node.io.filereader.FileAnalyzer;
import org.knime.base.node.io.filereader.FileReaderNodeSettings;
import org.knime.base.node.io.filereader.FileTable;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.NodeLogger;
import org.knime.core.util.tokenizer.SettingsStatus;

public class IO {

	private static FileWriterSettings getWriterSettings(){
		FileWriterSettings writerSettings     = new FileWriterSettings();
		
		// writer settings
		writerSettings.setWriteColumnHeader(true);
		writerSettings.setWriteRowID(true);
		writerSettings.setQuoteBegin("\"");
		writerSettings.setQuoteEnd("\"");	
		
		return writerSettings;
	}
	
	private static FileReaderNodeSettings getReaderSettings(boolean rowHeader, boolean colHeader){
		FileReaderNodeSettings readerSettings = new FileReaderNodeSettings();
		
		// reader settings
		readerSettings.addDelimiterPattern(";", false, false, false);
		readerSettings.setDelimiterUserSet(true);
		readerSettings.addRowDelimiter("\n", true);
		readerSettings.addQuotePattern("\"", "\"");
		readerSettings.setQuoteUserSet(true);
		readerSettings.addSingleLineCommentPattern("#", false, false);
		readerSettings.setCommentUserSet(true);
		readerSettings.setFileHasColumnHeaders(colHeader);
		readerSettings.setFileHasColumnHeadersUserSet(true);
		readerSettings.setFileHasRowHeaders(rowHeader);
		readerSettings.setFileHasRowHeadersUserSet(true);
		readerSettings.setWhiteSpaceUserSet(true);
		
		return readerSettings;
	}
	
	
	public static BufferedDataTable readCSV(final ExecutionContext exec,String filepath, boolean rowHeader, boolean colHeader) throws IOException, CanceledExecutionException{
		URL url = new URL("file:"+filepath);
		
		// set reader
		FileReaderNodeSettings readerSettings = getReaderSettings(rowHeader, colHeader);
		readerSettings.setDataFileLocationAndUpdateTableName(url);
		readerSettings = FileAnalyzer.analyze(readerSettings, null);
		SettingsStatus status = readerSettings.getStatusOfSettings();
		if (status.getNumOfErrors() > 0) {
			throw new IllegalStateException(status.getErrorMessage(0));
		}
		
		FileTable fTable = new FileTable(readerSettings.createDataTableSpec(),readerSettings, exec.createSubExecutionContext(0.5));
		BufferedDataTable table = exec.createBufferedDataTable(fTable, exec.createSubProgress(0.5));
		return table;
	}

	public static BufferedDataTable[] readCSV(ExecutionContext exec, String[] files, NodeLogger logger, boolean rowHeader, boolean colHeader) {
		BufferedDataTable[] result = new BufferedDataTable[files.length];

		for(int i=0; i<files.length;i++){
			try {
				result[i] = IO.readCSV(exec, files[i], rowHeader, colHeader);
			} catch (IOException e){
				logger.error("readCSV: IOException.");
				e.printStackTrace();
			} catch (CanceledExecutionException e ){
				logger.error("readCSV: Execution canceled.");
				e.printStackTrace();
			}
		}
		return result;
	}

	
	
	/**
	 * Writes a BufferedDataTable to file
	 * @param inData data to write
	 * @param file file to write to
	 * @param exec execution environment
	 * @throws CanceledExecutionException 
	 */
	public static void writeAsCSV(final BufferedDataTable inData, File file, ExecutionContext exec, NodeLogger logger) throws CanceledExecutionException {
		FileWriterSettings writerSettings = getWriterSettings();
		File parentDir = file.getParentFile();

		// create directory
		if (!parentDir.exists()) {
			if (!parentDir.mkdirs()) {
				logger.error("Unable to create directory for specified output file: " + parentDir.getAbsolutePath());
				return;
			}
			logger.info("Created directory for specified output file: " + parentDir.getAbsolutePath());
		}

		// create output stream
		OutputStream outStream = null;
		try {
			outStream = new FileOutputStream(file, false);
		} catch (FileNotFoundException e) {
			logger.error("Unable to create OutputStream");
			e.printStackTrace();
		}

		// table writer
		try {
			CSVWriter tableWriter = new CSVWriter(new OutputStreamWriter(outStream),writerSettings);
			tableWriter.write(inData, exec);
			tableWriter.close();

			if (tableWriter.hasWarningMessage()) {
				logger.warn(tableWriter.getLastWarningMessage());
			}

		} catch (IOException e) {
			logger.error("Writing to table '" + file.getAbsolutePath() + "' cancled");
			e.printStackTrace();
			throw(new CanceledExecutionException("Writing to table '" + file.getAbsolutePath() + "' cancled\n" + e.getMessage()));
		}
	}

	
	
	@Deprecated
	public static String findFile(String folder, String fileExtension, String curr_infile){
		//add, cv,lod,norm,plot,set
		//Alternative fuer alles ausser FileConverter und Info.file
		if(fileExtension.equals("add")&fileExtension.equals("cv")&fileExtension.equals("lod")&fileExtension.equals("norm")&fileExtension.equals("plot")&fileExtension.equals("set")){
			String base[] = curr_infile.split("csv");
			String newName = base[0]+fileExtension+".csv";
			System.out.println("#########################");
			System.out.println(folder+"/"+newName);
			System.out.println("#########################");
			return folder+"/"+newName;
		}
		//

		File[] files = new File(folder+"/").listFiles();
		fileExtension = fileExtension.replaceAll("\\+","\\\\+");
		fileExtension = fileExtension.replaceAll("\\*","\\\\*");
		for (File file : files) {
			if (!file.isDirectory()) {
				String name = file.getName();
				String [] namesplit = name.split(fileExtension);

				//System.out.println(fileExtension);
				if(namesplit[namesplit.length-1].equals(".csv")){
					System.out.println("File: " + folder+"/"+file.getName());
					return folder+"/"+name;
				}
			}
		}
		return "";
	}

	@Deprecated
	public static boolean deleteFile(String filepath){
		File f = new File(filepath);
		boolean suc = f.delete();
		return suc;
	}
}
