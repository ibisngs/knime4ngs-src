package de.helmholtz_muenchen.ibis.ngs.vqsr;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "VQSR" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class VQSRNodeDialog extends HTExecutorNodeDialog {

	
	/**
	 * Models
	 */
    
    private final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, "SNP");
    private final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0");
    private final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, "-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum");
    private final SettingsModelIntegerBounded m_GAUSSIANS = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_GAUSS,8,1,Integer.MAX_VALUE);
    private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

    private final SettingsModelString m_RESOURCES_STRING_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_1, "hapmap,known=false,training=true,truth=true,prior=15.0");
    private final SettingsModelString m_RESOURCES_STRING_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_2, "omni,known=false,training=true,truth=true,prior=12.0");
    private final SettingsModelString m_RESOURCES_STRING_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_3, "1000G,known=false,training=true,truth=false,prior=10.0");
    private final SettingsModelString m_RESOURCES_STRING_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_4, "dbsnp,known=true,training=false,truth=false,prior=2.0");
    private final SettingsModelString m_RESOURCES_STRING_5 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_5, "mills,known=false,training=true,truth=true,prior=12.0");
    
    private final SettingsModelString m_RESOURCES_FILE_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/hapmap_3.3.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_2, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_omni2.5.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_3, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_phase1.snps.high_confidence.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_4, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/dbsnp_138.hg19.vcf");
    private final SettingsModelString m_RESOURCES_FILE_5 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_5, "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf");

    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_1 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_2 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_2, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_3 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_3, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_4 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_4, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_5 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_5, false);

    private final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,99.9,1,100);
	
    private final SettingsModelOptionalString m_OPT_VAR_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_VAR_RECAL,"",false);
    private final SettingsModelOptionalString m_OPT_APPLY_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_APPLY_RECAL,"",false);
    /**
     * New pane for configuring the VQSR node.
     */
    protected VQSRNodeDialog() {
    	super();
    
    	createNewGroup("Path to GATK jar file");
    	addDialogComponent(new DialogComponentFileChooser(m_GATK, "gatk_vqsr", 0, ".jar"));
    	
    	createNewGroup("Reference Genome");
    	addDialogComponent(new DialogComponentFileChooser(m_REF_GENOME, "gatk_vqsr_ref_genome", 0, ".fa",".fasta",".txt"));
    	
    	createNewGroup("General Options");
    	addDialogComponent(new DialogComponentStringSelection(m_MODE, "Recalibration Mode", "SNP","INDEL"));
    	
    	createNewGroup("ApplyRecalibration");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_TS_FILTER, "TS Filter Level", 0.1,6));
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_APPLY_RECAL,"Optional flags"));
    	setHorizontalPlacement(false);
    	
    	createNewTab("VariantRecalibrator");
    	createNewGroup("General Options");
    	addDialogComponent(new DialogComponentString(m_TRANCHE, "Tranche Levels"));
    	addDialogComponent(new DialogComponentString(m_AN, "Annotation"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_GAUSSIANS, "Gaussians",1,4));
    	addDialogComponent(new DialogComponentNumber(m_NT, "Threads", 1,4));
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_VAR_RECAL,"Optional flags"));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Resources");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_1, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_1, ""));
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_1, "gatk_vqsr_resource_file1", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_2, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_2, ""));
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_2, "gatk_vqsr_resource_file2", 0, ".vcf"));  
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_3, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_3, ""));
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_3, "gatk_vqsr_resource_file3", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_4, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_4, ""));
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_4, "gatk_vqsr_resource_file4", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_5, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_5, ""));
    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_5, "gatk_vqsr_resource_file5", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	m_RESOURCES_BOOLEAN_1.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_RESOURCES_FILE_1.setEnabled(m_RESOURCES_BOOLEAN_1.getBooleanValue());
				m_RESOURCES_STRING_1.setEnabled(m_RESOURCES_BOOLEAN_1.getBooleanValue());
				
			}
		});
    	m_RESOURCES_BOOLEAN_2.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_RESOURCES_FILE_2.setEnabled(m_RESOURCES_BOOLEAN_2.getBooleanValue());
				m_RESOURCES_STRING_2.setEnabled(m_RESOURCES_BOOLEAN_2.getBooleanValue());	
			}
		});
    	m_RESOURCES_BOOLEAN_3.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_RESOURCES_FILE_3.setEnabled(m_RESOURCES_BOOLEAN_3.getBooleanValue());
				m_RESOURCES_STRING_3.setEnabled(m_RESOURCES_BOOLEAN_3.getBooleanValue());
			}
		});
    	m_RESOURCES_BOOLEAN_4.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_RESOURCES_FILE_4.setEnabled(m_RESOURCES_BOOLEAN_4.getBooleanValue());
				m_RESOURCES_STRING_4.setEnabled(m_RESOURCES_BOOLEAN_4.getBooleanValue());	
			}
		});
    	m_RESOURCES_BOOLEAN_5.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				m_RESOURCES_FILE_5.setEnabled(m_RESOURCES_BOOLEAN_5.getBooleanValue());
				m_RESOURCES_STRING_5.setEnabled(m_RESOURCES_BOOLEAN_5.getBooleanValue());	
			}
		});
    	
    	m_MODE.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				boolean isSNP = false;
				if(m_MODE.getStringValue().equals("SNP")) {
					isSNP = true;
				}
				m_RESOURCES_BOOLEAN_1.setBooleanValue(isSNP);
				m_RESOURCES_STRING_1.setEnabled(isSNP);
				m_RESOURCES_FILE_1.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_2.setBooleanValue(isSNP);
				m_RESOURCES_STRING_2.setEnabled(isSNP);
				m_RESOURCES_FILE_2.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_3.setBooleanValue(isSNP);
				m_RESOURCES_STRING_3.setEnabled(isSNP);
				m_RESOURCES_FILE_3.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_4.setBooleanValue(true);
				m_RESOURCES_STRING_4.setEnabled(true);
				m_RESOURCES_FILE_4.setEnabled(true);
				m_RESOURCES_BOOLEAN_5.setBooleanValue(!isSNP);
				m_RESOURCES_STRING_5.setEnabled(!isSNP);
				m_RESOURCES_FILE_5.setEnabled(!isSNP);
			}
    		
    	});
    }

	@Override
    protected void updatePrefs() {
    	if(usePrefPage.getBooleanValue()) {
    		String gatkPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("GenomeAnalysisTK.jar");
    		if(gatkPath != null && !gatkPath.equals("")) {
    			m_GATK.setStringValue(gatkPath);
    			m_GATK.setEnabled(false);
    		} else {
    			m_GATK.setEnabled(true);
    		}
    		String refGenome = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
    		if(refGenome != null && !refGenome.equals("")) {
    			m_REF_GENOME.setStringValue(refGenome);
    			m_REF_GENOME.setEnabled(false);
    		} else {
    			m_REF_GENOME.setEnabled(true);
    		}
    	} else {
    		m_GATK.setEnabled(true);
    		m_REF_GENOME.setEnabled(true);
    	}
    }
}

