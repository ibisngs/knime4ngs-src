package de.helmholtz_muenchen.ibis.ngs.snpeffgetdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.ngs.ShowOutput;

/**
 * This is the model implementation of SnpEffGetDB.
 * 
 *
 * @author 
 */
public class SnpEffGetDBNodeModel extends NodeModel {
    
	/**
	 * Input arguments
	 */
	public static final String CFGKEY_SNPEFF_FOLDER="snpeff_folder";
	public static final String CFGKEY_DATABASE="database";
	//public static final String CFGKEY_DATABASE_FOLDER="database_folder";
	
	//directory containing the snpEff.jar and scripts
	private final SettingsModelString m_snpeff_folder = new SettingsModelString(
			SnpEffGetDBNodeModel.CFGKEY_SNPEFF_FOLDER,"");
	//name of the database to download
	private final SettingsModelString m_database = new SettingsModelString(
			SnpEffGetDBNodeModel.CFGKEY_DATABASE,"");
	//directory where to the database is extracted
	//.config file will be modified!!!
	//private final SettingsModelString m_database_folder = new SettingsModelString(
	//		SnpEffGetDBNodeModel.CFGKEY_DATABASE_FOLDER,"");
	
	
    /**
     * Constructor for the node model.
     */
    protected SnpEffGetDBNodeModel() {
    
    	super(0,1);
    	
        // TODO: Specify the amount of input and output ports needed.

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	String snpEffDirectory = m_snpeff_folder.getStringValue();
    	String dbName = m_database.getStringValue();
    	//String dbFolder = m_database_folder.getStringValue();
    	
    	//TODO: change folder for log file or delete
    	/**Initialize logfile**/
    	String logfile = snpEffDirectory + "/logfile.txt";
    	ShowOutput.setLogFile(logfile);
    	StringBuffer logBuffer = new StringBuffer(50);
    	logBuffer.append(ShowOutput.getNodeStartTime("SnpEffGetDB"));
    	/**logfile initialized**/
    	
    	/*Modify .config file to set new database folder*/
    	String configFile = snpEffDirectory + "/snpEff.config";
    	String pattern = "^data_dir.*";
        Pattern p = Pattern.compile(pattern);
    	BufferedReader br = new BufferedReader(new FileReader(configFile));
        StringBuilder sb = new StringBuilder();
        
    	try {
            String line = br.readLine();
            while (line != null) {
                Matcher m = p.matcher(line);

                if (m.find( )) {
                	//sb.append("data_dir = " + dbFolder + "\n");
                	sb.append("data_dir = " + snpEffDirectory + "/databases\n");
                }
                else{
                	sb.append(line);
                }
                
                sb.append("\n");
                line = br.readLine();
            }

        } finally {
            br.close();
        }
    	//delete old config file and write new one
    	try{
    		File file = new File(configFile);
    		file.delete();
 
    	}catch(Exception e){
 
    		e.printStackTrace();
    	}
    	//write modified file
    	try {
    		FileWriter outFile = new FileWriter(configFile);
    		PrintWriter out = new PrintWriter(outFile);
    		out.print(sb.toString());
    		out.close();
    		} catch (IOException e){
        		e.printStackTrace();
        }
    	
    	/*Check if database does not exist yet
    	 * and download only if it does not exist*/
    	File file=new File(snpEffDirectory + "/databases/" + dbName);
    	if(!file.exists()) {
	        //Make and execute command
	    	//TODO: increase heap size for large databases?
	    	//Note: currently database is downloaded into the install directory.
	    	//		if this should be changed, add the -c option to specify the
	    	//		path to snpEff.config
	    	String command = "cd " + snpEffDirectory + "; java -jar snpEff.jar download -v " + dbName;
	    	//run the task
	      	ProcessBuilder b = new ProcessBuilder("/bin/sh", "-c", command);
	      	Process p_snpEff = b.start();
	    	p_snpEff.waitFor();  	
    	    	//TODO remove?
	    	logBuffer.append(ShowOutput.getLogEntry(p_snpEff, command));
	    	logBuffer.append(ShowOutput.getNodeEndTime());
	    	ShowOutput.writeLogFile(logBuffer);
    	}
    	else{
	    	logBuffer.append("Database " + dbName + " already exists\n");
	    	logBuffer.append(ShowOutput.getNodeEndTime());
	    	ShowOutput.writeLogFile(logBuffer);
    	}
    	
    	//make output table
    	DataColumnSpecCreator col1 = new DataColumnSpecCreator("snpEffDirectory", StringCell.TYPE);
        DataColumnSpecCreator col2 = new DataColumnSpecCreator("database", StringCell.TYPE);
        DataColumnSpec[] cols = new DataColumnSpec[]{col1.createSpec(),col2.createSpec()};
    	DataTableSpec table = new DataTableSpec(cols);
    	BufferedDataContainer cont = exec.createDataContainer(table);
    	StringCell cl1 = new StringCell(snpEffDirectory);
    	StringCell cl2 = new StringCell(dbName);
    	DataCell[] c = new DataCell[]{cl1,cl2};
    	DefaultRow r = new DefaultRow("Row0",c);
    	cont.addRowToTable(r);
    	cont.close();
    	BufferedDataTable out = cont.getTable();
    	    	
        return new BufferedDataTable[]{out};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        if(m_snpeff_folder.getStringValue().length() == 0){
        	throw new InvalidSettingsException("Specify path to snpEff!");
        }
        
        if(m_database.getStringValue().length() == 0){
        	throw new InvalidSettingsException("Specify a database to download!");
        }
        
        //if(m_database_folder.getStringValue().length() == 0){
        //	throw new InvalidSettingsException("Specify a database to download!");
        //}

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
        m_snpeff_folder.saveSettingsTo(settings);
    	m_database.saveSettingsTo(settings);
    	//m_database_folder.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        
        m_snpeff_folder.loadSettingsFrom(settings);
    	m_database.loadSettingsFrom(settings);
    	//m_database_folder.loadSettingsFrom(settings);
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
        m_snpeff_folder.validateSettings(settings);
    	m_database.validateSettings(settings);
    	//m_database_folder.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}

