package de.helmholtz_muenchen.ibis.ngs.filtersummary;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.knime.core.data.DataRow;
import org.knime.core.data.DataTableSpec;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;

import de.helmholtz_muenchen.ibis.ngs.filtersummary.vcf.LoFWorkflowResult;


/**
 * This is the model implementation of FilterSummary.
 * 
 *
 * @author 
 */
public class FilterSummaryNodeModel extends NodeModel {
    
    // the logger instance
    private static final NodeLogger logger = NodeLogger.getLogger(FilterSummaryNodeModel.class);
   
   // filter settings
   // minimum read coverage
   static final String CFGKEY_MIN_COVERAGE ="min_coverage";
   static final int DEF_MIN_COVERAGE=20;
   static final int MIN_MIN_COVERAGE=1;
   static final int MAX_MIN_COVERAGE=Integer.MAX_VALUE;
   private final SettingsModelIntegerBounded m_min_coverage = new SettingsModelIntegerBounded(CFGKEY_MIN_COVERAGE, DEF_MIN_COVERAGE, MIN_MIN_COVERAGE, MAX_MIN_COVERAGE);
   // minimum fraction of supporting reads
   static final String CFGKEY_SUPP_READS_FRAC="supp_reads";
   static final double DEF_SUPP_READS_FRAC=0.2;
   static final double MIN_SUPP_READS_FRAC=0;
   static final double MAX_SUPP_READS_FRAC=1;
   private final SettingsModelDoubleBounded m_supp_reads_frac = new SettingsModelDoubleBounded(CFGKEY_SUPP_READS_FRAC, DEF_SUPP_READS_FRAC, MIN_SUPP_READS_FRAC, MAX_SUPP_READS_FRAC);
   
   //summary settings
   // OMIM IDs for genes
   static final String CFGKEY_ANN_OMIM="ann_omim";
   static final boolean DEF_ANN_OMIM=false;
   private final SettingsModelBoolean m_ann_omim= new SettingsModelBoolean(CFGKEY_ANN_OMIM, DEF_ANN_OMIM);
   // gene info file
   static final String CFGKEY_GENE_INFO_FILE="gene_info_file";
   static final String DEF_GENE_INFO_FILE="";
   private final SettingsModelString m_gene_info_file = new SettingsModelString(CFGKEY_GENE_INFO_FILE, DEF_GENE_INFO_FILE);
   // pathways for genes
   static final String CFGKEY_ANN_PATHWAY="ann_pathway";
   static final boolean DEF_ANN_PATHWAY=false;
   private final SettingsModelBoolean m_ann_pathways = new SettingsModelBoolean(CFGKEY_ANN_PATHWAY, DEF_ANN_PATHWAY);
   // Wikipathway file
   static final String CFGKEY_WIKIPATHWAY_FILE="wikipathway_file";
   static final String DEF_WIKIPATHWAY_FILE="";
   private final SettingsModelString m_wikipathway_file = new SettingsModelString(CFGKEY_WIKIPATHWAY_FILE, DEF_WIKIPATHWAY_FILE);
   // number of interactions for genes
   static final String CFGKEY_ANN_NUM_PPI="ann_num_ppi";
   static final boolean DEF_ANN_NUM_PPI =false;
   private final SettingsModelBoolean m_ann_num_ppi = new SettingsModelBoolean(CFGKEY_ANN_NUM_PPI, DEF_ANN_NUM_PPI);
   // Reactome file
   static final String CFGKEY_REACTOME_FILE="reactome_file";
   static final String DEF_REACTOME_FILE="";
   private final SettingsModelString m_reactome_file = new SettingsModelString(CFGKEY_REACTOME_FILE, DEF_REACTOME_FILE);
   
   // statistics summary
   static final String CFGKEY_WIRTE_STATS="write_stats";
   static final boolean DEF_WRITE_STATS=true;
   private final SettingsModelBoolean m_write_stats = new SettingsModelBoolean(CFGKEY_WIRTE_STATS, DEF_WRITE_STATS);
   //output original vcf sample and info colum
   static final String CFGKEY_ORG_VCF_COL = "org_vcf_col";
   static final boolean DEF_ORG_VCF_COL=false;
   private final SettingsModelBoolean m_org_vcf_col= new SettingsModelBoolean(CFGKEY_ORG_VCF_COL, DEF_ORG_VCF_COL);
   
   private int posGATKSnp=-1;
   private int tableGATKSnp=-1;
   private int posGATKIndel=-1;
   private int tableGATKIndel=-1;
   private int posPindelD=-1;
   private int tablePindelD=-1;
   private int posPindelSI=-1;
   private int tablePindelSI=-1;
   
    /**
     * Constructor for the node model.
     */
    protected FilterSummaryNodeModel() {
    
        // 2 incoming, no outcoming    	
        super(2,0);
        
        m_gene_info_file.setEnabled(false);
        m_wikipathway_file.setEnabled(false);
        m_reactome_file.setEnabled(false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
        
        DataRow [] r = {inData[0].iterator().next(), inData[1].iterator().next()};
        
        String pathgatksnps=r[tableGATKSnp].getCell(posGATKSnp).toString();
        String pathgatkindel=r[tableGATKIndel].getCell(posGATKIndel).toString();
        String pathpindeld=r[tablePindelD].getCell(posPindelD).toString();
        String pathpindelsi=r[tablePindelSI].getCell(posPindelSI).toString();
      
        // check if path is null
        if(pathgatksnps.equals("")){
        	throw new Exception("No GATK snp vcf file available, something went wrong with the previous node!");
        }
        if(pathgatkindel.equals("")){
        	throw new Exception("No GATK indel vcf file available, something went wrong with the previous node!");
        }
        if(pathpindeld.equals("")){
        	throw new Exception("No Pindel deletion vcf file available, something went wrong with the previous node!");
        }
        if(pathpindelsi.equals("")){
        	throw new Exception("No Pindel insertion vcf file available, something went wrong with the previous node!");
        }
        
        // check path to vcf file
        if(!Files.exists(Paths.get(pathgatksnps))){
        	throw new Exception("Path to input GATK snp vcf file: "+pathgatksnps+" does not exist");
        }
        if(!Files.exists(Paths.get(pathgatkindel))){
        	throw new Exception("Path to input GATK indel vcf file: "+pathgatkindel+" does not exist");
        }
        if(!Files.exists(Paths.get(pathpindeld))){
        	throw new Exception("Path to input GATK snp vcf file: "+pathpindeld+" does not exist");
        }
        if(!Files.exists(Paths.get(pathpindelsi))){
        	throw new Exception("Path to input GATK snp vcf file: "+pathpindelsi+" does not exist");
        }
        
      	String wikipath=null;
      	String reactome=null;
      	String geneinfo=null;
      	
      	if(m_ann_omim.getBooleanValue()){
      		geneinfo=m_gene_info_file.getStringValue();
      		if(geneinfo.equals("")){
      			throw new Exception("Missing gene_info file: You have to configure the node before executing it!");
      		}
      	    if(!Files.exists(Paths.get(geneinfo))){
      	    	throw new Exception("Gene_info file: "+geneinfo+" does not exist");
      	    }
      	}
      	
      	if(m_ann_pathways.getBooleanValue()){
      		wikipath=m_wikipathway_file.getStringValue();
      		if(wikipath.equals("")){
      			throw new Exception("Missing wikipathway file: You have to configure the node before executing it!");
      		}
      	    if(!Files.exists(Paths.get(wikipath))){
      	    	throw new Exception("Wikipathway file: "+wikipath+" does not exist");
      	    }
      	}
      	
      	if(m_ann_num_ppi.getBooleanValue()){
      		reactome=m_reactome_file.getStringValue();
      		if(reactome.equals("")){
      			throw new Exception("Missing Reactome file: You have to configure the node before executing it!");
      		}
      	    if(!Files.exists(Paths.get(reactome))){
      	    	throw new Exception("Reactome file: "+reactome+" does not exist");
      	    }
      	}
      	
    	
      	// create path to for outputfiles
      	String outbase="";
      	Path p = Paths.get(Paths.get(pathgatksnps).getParent().toString(), "filtersummary");
		//file does not exist
		if(!Files.exists(Paths.get(p.toString()+".all.csv"))){
			outbase=p.toString();
		}
		//file exists
		else{
			
			//add integer to file name and increment integer until the file does not exist
			int n=1;
			p= Paths.get(Paths.get(pathgatksnps).getParent().toString(), "filtersummary"+n);
			while(Files.exists(Paths.get(p.toString()+".all.csv"))){
				n++;

				if(n==10000){
					throw new Exception("Oops, I'm so sorry, something went wrong, there are too many files");
				}
				p= Paths.get(Paths.get(pathgatksnps).getParent().toString(), "filtersummary"+n);
			}
			outbase=p.toString();
		}
        outbase+=".";
        
        // run 
		LoFWorkflowResult lwr = new LoFWorkflowResult(
				pathpindeld,
				pathpindelsi, 
				pathgatksnps, 
				pathgatkindel, 
				outbase,
				geneinfo,
				wikipath,
				reactome);
		
		lwr.cf.setFilterOptions(m_min_coverage.getIntValue(), m_supp_reads_frac.getDoubleValue());
		
		int original_cols = 1;
		if(m_org_vcf_col.getBooleanValue()){
			original_cols=2;
		}
		
		lwr.runAnalysis(m_write_stats.getBooleanValue(), original_cols, false);

        return null;
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
        
    	/*requires input of 4 vcf files
    	 * Path2VCFsnpFile -> GATK
    	 * Path2VCFindelFile -> GATK
    	 * Path2VCFinsertionsFile -> Pindel
    	 * Path2VCFdeletionsFile -> Pindel
    	 */
    	
    	DataTableSpec inport1 = inSpecs[0];
    	DataTableSpec inport2 = inSpecs[1];
    	
    	if(!inport1.containsName("Path2VCFsnpFile") && !inport2.containsName("Path2VCFsnpFile")){
    		throw new InvalidSettingsException("Missing GATK snps file from previous node!");
    	}
    	if(!inport1.containsName("Path2VCFindelFile") && !inport2.containsName("Path2VCFindelFile")){
    		throw new InvalidSettingsException("Missing GATK indel file from previous node!");
    	}
    	if(!inport1.containsName("Path2VCFinsertionsFile") && !inport2.containsName("Path2VCFinsertionsFile")){
    		throw new InvalidSettingsException("Missing Pindel insertion file from previous node!");
    	}
    	if(!inport1.containsName("Path2VCFdeletionsFile") && !inport2.containsName("Path2VCFdeletionsFile")){
    		throw new InvalidSettingsException("Missing Pindel deletion file from previous node!");
    	}
    	
    	String [] names1 = inport1.getColumnNames();
    	String [] names2 = inport2.getColumnNames();
    	
    	for(int i=0; i<names1.length; i++){
    		if(names1[i].equals("Path2VCFsnpFile")){
    			posGATKSnp=i;
    			tableGATKSnp=0;
    		}
    		if(names1[i].equals("Path2VCFindelFile")){
    			posGATKIndel=i;
    			tableGATKIndel=0;
    		}
    		if(names1[i].equals("Path2VCFinsertionsFile")){
    			posPindelSI=i;
    			tablePindelSI=0;
    		}
    		if(names1[i].equals("Path2VCFdeletionsFile")){
    			posPindelD=i;
    			tablePindelD=0;
			}
    	}
    	
    	for(int i=0; i<names2.length; i++){
    		if(names2[i].equals("Path2VCFsnpFile")){
    			posGATKSnp=i;
    			tableGATKSnp=1;
    		}
    		if(names2[i].equals("Path2VCFindelFile")){
    			posGATKIndel=i;
    			tableGATKIndel=1;
    		}
    		if(names2[i].equals("Path2VCFinsertionsFile")){
    			posPindelSI=i;
    			tablePindelSI=1;
    		}
    		if(names2[i].equals("Path2VCFdeletionsFile")){
    			posPindelD=i;
    			tablePindelD=1;
			}
    	}

        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
        
        m_min_coverage.saveSettingsTo(settings);
        m_supp_reads_frac.saveSettingsTo(settings);
        m_write_stats.saveSettingsTo(settings);
        m_org_vcf_col.saveSettingsTo(settings);
        m_ann_omim.saveSettingsTo(settings);
        m_gene_info_file.saveSettingsTo(settings);
        m_ann_pathways.saveSettingsTo(settings);
        m_wikipathway_file.saveSettingsTo(settings);
        m_ann_num_ppi.saveSettingsTo(settings);
        m_reactome_file.saveSettingsTo(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {

    	m_min_coverage.loadSettingsFrom(settings);
    	m_supp_reads_frac.loadSettingsFrom(settings);
    	m_write_stats.loadSettingsFrom(settings);
    	m_org_vcf_col.loadSettingsFrom(settings);
    	m_ann_omim.loadSettingsFrom(settings);
    	m_gene_info_file.loadSettingsFrom(settings);
    	m_ann_pathways.loadSettingsFrom(settings);
    	m_wikipathway_file.loadSettingsFrom(settings);
    	m_ann_num_ppi.loadSettingsFrom(settings);
    	m_reactome_file.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
    
    	m_min_coverage.validateSettings(settings);
    	m_supp_reads_frac.validateSettings(settings);
    	m_write_stats.validateSettings(settings);
    	m_org_vcf_col.validateSettings(settings);
    	m_ann_omim.validateSettings(settings);
    	m_gene_info_file.validateSettings(settings);
    	m_ann_pathways.validateSettings(settings);
    	m_wikipathway_file.validateSettings(settings);
    	m_ann_num_ppi.validateSettings(settings);
    	m_reactome_file.validateSettings(settings);
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

