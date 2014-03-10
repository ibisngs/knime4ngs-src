package de.helmholtz_muenchen.ibis.ngs.vat;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.DefaultRow;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.utils.lofs.PathProcessor;


/**
 * This is the model implementation of VAT.
 * 
 *
 * @author 
 */
public class VATNodeModel extends NodeModel {
    
    // the logger instance
    protected static final NodeLogger logger = NodeLogger.getLogger(VATNodeModel.class);
        
    /** the settings key which is used to retrieve and 
        store the settings (from the dialog or from a settings file)    
       (package visibility to be usable from the dialog). */
    
    // path to folder containing VAT executables: snpMapper, indelMapper
    static final String CFGKEY_VAT_FOLDER="vat_folder";
    static final String DEF_VAT_FOLDER="";
	private final SettingsModelString m_VAT_folder = new SettingsModelString(CFGKEY_VAT_FOLDER, DEF_VAT_FOLDER);
	// path to files gene intervals -> genomic positions of exons
	static final String CFGKEY_INTERVALS="intervals";
	static final String DEF_INTERVALS="";
	private final SettingsModelString m_intervals = new SettingsModelString(CFGKEY_INTERVALS, DEF_INTERVALS);
	// path to file with transcript sequences
	static final String CFGKEY_TRANSCRIPTS="transcripts";
	static final String DEF_TRANSCRIPTS="";
	private final SettingsModelString m_transcripts = new SettingsModelString(CFGKEY_TRANSCRIPTS, DEF_TRANSCRIPTS);
    
	// if files are available from previous nodes and their positions
	private boolean insertions;		//pindel_SI
	private boolean deletions;		//pindel_D
	private boolean snps;			//gatk snps
	private boolean indels;			//gatk indels
	
	private int posInsertions;
	private int posDeletions;
	private int posSnps;
	private int posIndels;
	
	

    /**
     * Constructor for the node model.
     */
    protected VATNodeModel() {
    
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {

    	// check data from previous node
    	
        //retrieve information from table
        DataRow r=inData[0].iterator().next();
        
        // check vcf files
        String deletionfile="";
        String insertionfile="";
        String snpfile="";
        String indelfile="";
        
        if(deletions){
        	
	        deletionfile=r.getCell(posDeletions).toString();
	        
	        // check if path is null
	        if(deletionfile.equals("")){
	        	throw new Exception("No vcf file with deletions available, something went wrong with the previous node!");
	        }
	        
	        // check path to file
	        if(!Files.exists(Paths.get(deletionfile))){
	        	throw new Exception("Path to vcf file with deletions: "+deletionfile+" does not exist!");
	        }
	        
	        if(!PathProcessor.getExt(deletionfile).equals("vcf")){
	        	throw new Exception("Input file: "+deletionfile+" is not in vcf format!");
	        }
        }
        
        if(insertions){
        	
	        insertionfile=r.getCell(posInsertions).toString();
	        
	        // check if path is null
	        if(insertionfile.equals("")){
	        	throw new Exception("No vcf file with insertions available, something went wrong with the previous node!");
	        }
	        
	        // check path to file
	        if(!Files.exists(Paths.get(insertionfile))){
	        	throw new Exception("Path to vcf file with insertions: "+insertionfile+" does not exist!");
	        }
	        
	        if(!PathProcessor.getExt(insertionfile).equals("vcf")){
	        	throw new Exception("Input file: "+insertionfile+" is not in vcf format!");
	        }
        }
        
        if(snps){
        	
	        snpfile=r.getCell(posSnps).toString();
	        
	        // check if path is null
	        if(snpfile.equals("")){
	        	throw new Exception("No vcf file with SNPs available, something went wrong with the previous node!");
	        }
	        
	        // check path to file
	        if(!Files.exists(Paths.get(snpfile))){
	        	throw new Exception("Path to vcf file with SNPs: "+snpfile+" does not exist!");
	        }
	        
	        if(!PathProcessor.getExt(snpfile).equals("vcf")){
	        	throw new Exception("Input file: "+snpfile+" is not in vcf format!");
	        }
        }
        
        if(indels){
        	
	        indelfile=r.getCell(posIndels).toString();
	        
	        // check if path is null
	        if(indelfile.equals("")){
	        	throw new Exception("No vcf file with Indels available, something went wrong with the previous node!");
	        }
	        
	        // check path to file
	        if(!Files.exists(Paths.get(indelfile))){
	        	throw new Exception("Path to vcf file with Indels: "+indelfile+" does not exist!");
	        }
	        
	        if(!PathProcessor.getExt(indelfile).equals("vcf")){
	        	throw new Exception("Input file: "+indelfile+" is not in vcf format!");
	        }
        }
        
        // check node settings
        
        //gene intervals
        String intervalfile=m_intervals.getStringValue();
        
        // check if path is null
        if(intervalfile.equals("")){
        	throw new Exception("Missing path to gene interval file: You have to configure the node before executing it!");
        }
        
        // check path to file
        if(!Files.exists(Paths.get(intervalfile))){
        	throw new Exception("Path to gene interval file: "+intervalfile+" does not exist!");
        }
        
        // fasta file with transcripts
        String transcriptfile=m_transcripts.getStringValue();
        
        // check if path is null
        if(transcriptfile.equals("")){
        	throw new Exception("Missing path to transcript file: You have to configure the node before executing it!");
        }
        
        // check path to file
        if(!Files.exists(Paths.get(transcriptfile))){
        	throw new Exception("Path to transcript file: "+transcriptfile+" does not exist!");
        }
        
        //VAT executables
        String vat=m_VAT_folder.getStringValue();
        
        // check if path is null
        if(vat.equals("")){
        	throw new Exception("Missing path to VAT folder: You have to configure the node before executing it!");
        }
        
        // check path to file
        if(!Files.exists(Paths.get(vat))){
        	throw new Exception("Path to VAT folder: "+vat+" does not exist!");
        }
        
        //check snpMapper
        String snpMapper="";
        
        if(snps){
        	if(!Files.exists(Paths.get(vat, "snpMapper"))){
        		throw new Exception("VAT folder does not contain snpMapper executable: "+Paths.get(vat, "snpMapper"));
        	}
        	else{
        		snpMapper=Paths.get(vat, "snpMapper").toString();
        	}
        }
        
        //check indelMapper
        String indelMapper="";
        
        if(deletions || insertions || indels){
        	if(!Files.exists(Paths.get(vat, "indelMapper"))){
        		throw new Exception("VAT folder does not contain indelMapper executable: "+Paths.get(vat, "indelMapper"));
        	}
        	else{
        		indelMapper=Paths.get(vat, "indelMapper").toString();
        	}
        }
        
        String deletionsout="";
        if(deletions){
        	
        	String deletionbase = PathProcessor.getBase(deletionfile);
        	deletionsout=PathProcessor.createOutputFile(deletionbase, "vcf", "vat");
        	
        	RunVAT.IndelMapper(exec, indelMapper, intervalfile, transcriptfile, deletionfile, deletionsout);
        }
        
        String insertionsout="";
        if(insertions){
        	
        	String insertionbase = PathProcessor.getBase(insertionfile);
        	insertionsout=PathProcessor.createOutputFile(insertionbase, "vcf", "vat");
        	
        	RunVAT.IndelMapper(exec, indelMapper, intervalfile, transcriptfile, insertionfile, insertionsout);
        }
        
        String indelsout="";
        if(indels){
        	
        	String indelbase = PathProcessor.getBase(indelfile);
        	indelsout=PathProcessor.createOutputFile(indelbase, "vcf", "vat");
        	
        	RunVAT.IndelMapper(exec, indelMapper, intervalfile, transcriptfile, indelfile, indelsout);
        }
        
        String snpsout="";
        if(snps){
        	
        	String snpbase = PathProcessor.getBase(snpfile);
        	snpsout=PathProcessor.createOutputFile(snpbase, "vcf", "vat");
        	
        	RunVAT.SNPMapper(exec, snpMapper, intervalfile, transcriptfile, snpfile, snpsout);
        }
        
        //create output table
        
        int count =0;
        if(deletions){
        	count++;
        }
        if(insertions){
        	count++;
        }
        if(indels){
        	count++;
        }
        if(snps){
        	count++;
        }
        
        // should not be necessary
        if(count == 0){
        	return inData;
        }
        
        DataColumnSpec [] colspec = new DataColumnSpec[count];
        
        int pos=0;
    	if(deletions){
	    	colspec[pos++]=new DataColumnSpecCreator("Path2VCFdeletionsFile", StringCell.TYPE).createSpec();    		
    	}
    	if(insertions){
    		colspec[pos++]=new DataColumnSpecCreator("Path2VCFinsertionsFile", StringCell.TYPE).createSpec();
    	}
    	if(snps){
    		colspec[pos++]=new DataColumnSpecCreator("Path2VCFsnpFile", StringCell.TYPE).createSpec();
    	}
    	if(indels){
    		colspec[pos++]=new DataColumnSpecCreator("Path2VCFindelFile", StringCell.TYPE).createSpec();
    	}
    	
    	//create table
	    DataTableSpec outspec=new DataTableSpec(colspec);
	    BufferedDataContainer c = exec.createDataContainer(outspec);
	    
	    // fill string cells
	    pos=0;
	    DataCell [] row = new DataCell [count];
	    if(deletions){
	    	row[pos++]=new StringCell(deletionsout);
	    }
	    if(insertions){
	    	row[pos++]=new StringCell(insertionsout);
	    }
	    if(snps){
	    	row[pos++]=new StringCell(snpsout);
	    }
	    if(indels){
	    	row[pos++]=new StringCell(indelsout);
	    }

	    //create row and add it to the container
	    c.addRowToTable(new DefaultRow("row0", row));
	    
	    //create final table
	    c.close();
	    BufferedDataTable out=c.getTable();

        return new BufferedDataTable [] {out};

    }

    /**
     * {@inheritDoc}
     */ 
    @Override
    protected void reset() {
        // Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
    	if(!inSpecs[0].containsName("Path2VCFdeletionsFile")
    			&& !inSpecs[0].containsName("Path2VCFinsertionsFile") 
    			&& !inSpecs[0].containsName("Path2VCFsnpFile") 
    			&& !inSpecs[0].containsName("Path2VCFindelFile")){
    		
    		throw new InvalidSettingsException("Previous node is incompatibe with this node: Missing vcf file");
    		
    	}
    	
    	String [] cols = inSpecs[0].getColumnNames();
    	
    	for (int i=0; i<cols.length; i++){
    		
    		// deletions from pindel
    		if(cols[i].equals("Path2VCFdeletionsFile")){
    			deletions=true;
    			posDeletions=i;
    		}
    		// insertions from pindel
    		if(cols[i].equals("Path2VCFinsertionsFile")){
    			insertions=true;
    			posInsertions=i;
    		}
    		// snps from gatk
    		if(cols[i].equals("Path2VCFsnpFile")){
    			snps=true;
    			posSnps=i;
    		}
    		//indels from gatk
    		if(cols[i].equals("Path2VCFindelFile")){
    			indels=true;
    			posIndels=i;
    		}   		
    		
    	}

        return new DataTableSpec[]{null};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

       m_VAT_folder.saveSettingsTo(settings);
       m_intervals.saveSettingsTo(settings);
       m_transcripts.saveSettingsTo(settings); 

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
    	
    	m_VAT_folder.loadSettingsFrom(settings);
    	m_intervals.loadSettingsFrom(settings);
    	m_transcripts.loadSettingsFrom(settings);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    	
    	m_VAT_folder.validateSettings(settings);
    	m_intervals.validateSettings(settings);
    	m_transcripts.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }

}

