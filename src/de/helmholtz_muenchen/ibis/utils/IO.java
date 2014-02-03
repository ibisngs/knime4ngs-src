package de.helmholtz_muenchen.ibis.utils;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

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

	/**
	 * Read a single csv file as Buffered Table
	 * @param exec
	 * @param filepath
	 * @param rowHeader
	 * @param colHeader
	 * @return
	 * @throws IOException
	 * @throws CanceledExecutionException
	 */
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

	/**
	 * Reads csv files as Buffered Tables
	 * @param exec
	 * @param files
	 * @param logger
	 * @param rowHeader
	 * @param colHeader
	 * @return
	 */
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

	//////////////////////////////////////////////////////////////////////////
	// GET SCRIPT PATH FOR EXECUTING E.G. R-SCRIPTS
	//////////////////////////////////////////////////////////////////////////
	private static String SCRIPT_PATH;
	private static void initScriptPath() {
		if(SCRIPT_PATH != null) return;
		// LOCALLY OR FROM JAR?
		String path = IO.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		if(path.endsWith("jar")){
			SCRIPT_PATH = "";
			try {
				String tmpFolder = createTempDirectory("knime_libs_", "").getCanonicalPath();
				IO.unzipJar(path, new File(tmpFolder));
				SCRIPT_PATH = tmpFolder + File.separatorChar;
			} catch (IOException e) {
				System.err.println("Couldn't unpack jar to tmp directory");
				e.printStackTrace();
			}
		}else{
			// EXECUTED LOCALLY
			SCRIPT_PATH = path + File.separatorChar;
		}
	}

	public static String getScriptPath() {
		if(SCRIPT_PATH == null){
			initScriptPath();
		}
		return(SCRIPT_PATH);
	}

	/**
	 * Unzip jar file to directory
	 * @param jarFile path to jar file
	 * @param destDir path to target directory
	 * @throws IOException
	 */
	private static void unzipJar(String jarFile, File destDir) throws IOException{
		if (!destDir.exists()) {
			System.out.println("JONAS: creating directory: " + destDir.getCanonicalPath());
			boolean result = destDir.mkdir();  

			if(!result) {    
				throw(new IOException("Couldn't create temporary directory to extract scripts"));
			}
		}

		JarFile jar = new JarFile(jarFile);
		Enumeration<JarEntry> files = jar.entries();
		while (files.hasMoreElements()) {
			JarEntry file = files.nextElement();
			File f = new File(destDir.getCanonicalPath() + File.separator + file.getName());
			if (file.isDirectory()) { // if its a directory, create it
				f.mkdir();
				continue;
			}
			InputStream is = jar.getInputStream(file); // get the input stream
			FileOutputStream fos = new FileOutputStream(f);
			while (is.available() > 0) {  // write contents of 'is' to 'fos'
				fos.write(is.read());
			}
			fos.close();
			is.close();
		}
	}

	/**
	 * Create a temporary directory in system's standard tmp dir
	 * @param prefix prefix for directory name
	 * @param suffix suffix for directory name
	 * @return File object of the created tmp dir
	 * @throws IOException
	 */
	public static File createTempDirectory(String prefix, String suffix) throws IOException {
		final File temp = File.createTempFile(prefix, suffix);

		if(!(temp.delete())){
			throw new IOException("Could not delete temp file: " + temp.getAbsolutePath());
		}
		if(!(temp.mkdirs())){
			throw new IOException("Could not create temp directory: " + temp.getAbsolutePath());
		}
		return (temp);
	}

	/**
	 * returns the last x lines of the given file as String
	 * @param src source file
	 * @param maxLines number of lines to return
	 * @return last maxLines Lines of the src file
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static String tail(File src, int maxLines) throws FileNotFoundException, IOException {
		BufferedReader reader = new BufferedReader(new FileReader(src));
		String[] lines = new String[maxLines];
		int lastNdx = 0;
		String out = "";
		for (String line=reader.readLine(); line != null; line=reader.readLine()) {
			if (lastNdx == lines.length) {
				lastNdx = 0;
			}
			lines[lastNdx++] = line;
		}

		for (int ndx=lastNdx; ndx != lastNdx-1; ndx++) {
			if (ndx == lines.length) {
				ndx = 0;
			}
			out+=lines[ndx]+"\n";
		}

		return out;
	}
}
