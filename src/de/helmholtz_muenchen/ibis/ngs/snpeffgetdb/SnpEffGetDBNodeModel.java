package de.helmholtz_muenchen.ibis.ngs.snpeffgetdb;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.SuccessfulRunChecker;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeModel;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCell;
import de.helmholtz_muenchen.ibis.utils.datatypes.file.FileCellFactory;
import de.helmholtz_muenchen.ibis.utils.threads.Executor;

/**
 * This is the model implementation of SnpEffGetDB.
 * 
 *
 * @author 
 */
public class SnpEffGetDBNodeModel extends HTExecutorNodeModel {
    
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
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(SnpEffGetDBNodeModel.class);
	
	//The Output Col Names
	public static final String OUT_COL1 = "snpEffDirectory";
	public static final String OUT_COL2 = "database";
	
	
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
    		ArrayList<String> command = new ArrayList<String>();
    		command.add("java");
    		command.add("-jar snpEff.jar download");
        	command.add("-v "+dbName);
        	
        	/**Execute**/
        	String lockFile = file.getAbsolutePath() + SuccessfulRunChecker.LOCK_ENDING;
        	super.executeCommand(new String[]{StringUtils.join(command, " ")}, exec, new File(lockFile));
        	
        	Executor.executeCommand(new String[]{StringUtils.join(command, " ")},exec,LOGGER);

    	}
    	else{
	    	System.out.println("Database " + dbName + " already exists\n");
    	}
    	
    	/**
    	 * Output
    	 */
    	BufferedDataContainer cont = exec.createDataContainer(
    			new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()}));
    	
    	FileCell[] c = new FileCell[]{
    			(FileCell) FileCellFactory.create(snpEffDirectory),
    			(FileCell) FileCellFactory.create(dbName)};
    	
    	cont.addRowToTable(new DefaultRow("Row0",c));
    	cont.close();
    	BufferedDataTable outTable = cont.getTable();
    	    	
        return new BufferedDataTable[]{outTable};
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

        return new DataTableSpec[]{new DataTableSpec(
    			new DataColumnSpec[]{
    					new DataColumnSpecCreator(OUT_COL1, FileCell.TYPE).createSpec(),
    					new DataColumnSpecCreator(OUT_COL2, FileCell.TYPE).createSpec()})};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
    	super.saveSettingsTo(settings);
    	    	
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
        
    	super.loadValidatedSettingsFrom(settings);
    	
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
    	
    	super.validateSettings(settings);
    	
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

