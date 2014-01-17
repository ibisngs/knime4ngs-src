package de.helmholtz_muenchen.ibis.utils;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.port.PortType;

public class Utils {

	/**
	 * 
	 * @param nrDataPorts
	 * @param optionalPortsIds
	 * @return
	 */
	@Deprecated
    public static PortType[] createOptionalPorts(final int nrDataPorts, final int... optionalPortsIds){
        PortType[] portTypes = new PortType[nrDataPorts];
        Arrays.fill(portTypes, BufferedDataTable.TYPE);        

        if (optionalPortsIds.length > 0) {
            for (int portId : optionalPortsIds) {
                if ((portId - 1) < nrDataPorts) {
                    portTypes[portId - 1] = new PortType(BufferedDataTable.class, true);;
                }
            }
        }
        return portTypes;
    } 
    
	/**
	 * executes command
	 * @param com command to execute
	 * @param logBuffer StringBuffer to append stdout and stderr of command
	 * @throws Exception
	 */
	public static void executeCommand(String com, StringBuffer logBuffer) throws Exception{
	  	System.out.println("EXECUTING: >" +com);
	 	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
		Process p1 = b.start();
		p1.waitFor();
		logBuffer.append("STDERR:\n");
		logBuffer.append(ShowOutput.getLogEntryStdErr(p1, com));
		if(p1.exitValue()!=0){
			String logtail = Utils.tail(logBuffer, 50);
			throw new Exception("Encountered error while executing the command!\nChecking logfile for details...\n"+logtail);
		}
		logBuffer.append("\n\n\nSTDOUT:\n");
		logBuffer.append(ShowOutput.getLogEntryStdOut(p1, com));
	}
	
	/**
	 * executes command
	 * @param com command to execute
	 * @param logBufferstdOut StringBuffer to append stdout
	 * @param logBufferstdErrStringBuffer to append stderr
	 * @throws Exception
	 */
	public static void executeCommand(String com,StringBuffer logBufferstdOut, StringBuffer logBufferstdErr, final NodeLogger logger) throws Exception{
	  	System.out.println(com);

	  	// Build Process
	  	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", com);
	  	Process p1;
		p1 = b.start();
		p1.waitFor();
		logBufferstdOut.append(ShowOutput.getLogEntryStdOut(p1, com));
		logBufferstdErr.append(ShowOutput.getLogEntryStdErr(p1, com));

		if(p1.exitValue()!=0){
			//Get last lines of logfile
			String logtail = logBufferstdErr.substring(Math.max(0,logBufferstdErr.length()-800), logBufferstdErr.length());
			logger.error("Error while executing command >"+com+"\nSTDERR:\n" + logBufferstdErr + "STDOUT:\n" + logBufferstdOut);
			throw new Exception("Encountered error while executing the command!\nChecking logfile for details...\n"+logtail);
		}
	}
	/**
	 * Returns last X lines of the given StringBuffer
	 * @param buf StringBuffer
	 * @param lines number of lines to return
	 * @return last X lines of the given StringBuffer
	 */
	@Deprecated
	public static String tail(StringBuffer buf, int lines){
		int idx = buf.length();
		while(lines>0 && idx>=0){
			idx = buf.lastIndexOf(System.getProperty("line.separator"), idx);
		}
		return buf.substring(idx+1);
	}
	
	
	/**
	 * returns the last x lines of the given file as String
	 * @param src source file
	 * @param maxLines number of lines to return
	 * @return last maxLines Lines of the src file
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	@Deprecated
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
	
	/**
	 * Create a temporary directory in system's standard tmp dir
	 * @param prefix prefix for directory name
	 * @param suffix suffix for directory name
	 * @return File object of the created tmp dir
	 * @throws IOException
	 */
	@Deprecated
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
	 * Returns timestamp of current system time
	 * @return time as string (format: "yyyy-MM-dd HH:mm:ss")
	 */
	public static String getCurrentTimeStamp() {
		SimpleDateFormat sdfDate = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");//dd/MM/yyyy
		Date now = new Date();
		String strDate = sdfDate.format(now);
		return strDate;
	}
	
	
	//////////////////////////////////////////////////////////////////////////
	// GET SCRIPT PATH FOR EXECUTING E.G. R-SCRIPTS
	//////////////////////////////////////////////////////////////////////////
	private static final String TMP_FOLDER = "/tmp/knime_libs/";
	private static String SCRIPT_PATH;
	private static void initScriptPath(){
		if(SCRIPT_PATH != null) return;
		// LOCALLY OR FROM JAR?
		String path = Utils.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		if(path.endsWith("jar")){
			// EXECUTED FROM JAR
			
			StringBuffer log = new StringBuffer();
			String com = "unzip -n "+ path + " -d " + TMP_FOLDER;
			try {
				Utils.executeCommand(com, log);
			} catch (Exception e) {
				e.printStackTrace();
			}
			SCRIPT_PATH = TMP_FOLDER + File.separatorChar;
		}else{
			// EXECUTED LOCALLY
			SCRIPT_PATH = path + File.separatorChar;
		}
	}
	
	public static String getScriptPath(){
		if(SCRIPT_PATH == null){
			initScriptPath();
		}
		return(SCRIPT_PATH);
		
	}
	
}
