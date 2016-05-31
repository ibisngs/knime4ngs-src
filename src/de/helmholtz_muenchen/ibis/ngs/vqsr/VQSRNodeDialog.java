package de.helmholtz_muenchen.ibis.ngs.vqsr;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
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
     * New pane for configuring the VQSR node.
     */
    protected VQSRNodeDialog() {}
    
    public void addToolDialogComponents() {
    	
    	/**
    	 * Models
    	 */
        
        final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
        final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
        final SettingsModelIntegerBounded m_XMX = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_XMX, 8, 1, Integer.MAX_VALUE);
        final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, VQSRNodeModel.DEFAULT_MODE);
        final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, VQSRNodeModel.DEFAULT_TRANCHES);
        final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, VQSRNodeModel.DEFAULT_SNP_AN);
        final SettingsModelIntegerBounded m_GAUSSIANS = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_GAUSS,8,1,Integer.MAX_VALUE);
        final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

        final SettingsModelString m_RESOURCES_STRING_HAPMAP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_HAPMAP, "hapmap,known=false,training=true,truth=true,prior=15.0");
        final SettingsModelString m_RESOURCES_STRING_OMNI = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_OMNI, "omni,known=false,training=true,truth=true,prior=12.0");
        final SettingsModelString m_RESOURCES_STRING_1000G = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_1000G, "1000G,known=false,training=true,truth=false,prior=10.0");
        final SettingsModelString m_RESOURCES_STRING_DBSNP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_DBSNP, "dbsnp,known=true,training=false,truth=false,prior=2.0");
        final SettingsModelString m_RESOURCES_STRING_MILLS = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCE_MILLS, "mills,known=false,training=true,truth=true,prior=12.0");
        
        final SettingsModelString m_RESOURCES_FILE_HAPMAP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_HAPMAP,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/hapmap_3.3.hg19.vcf");
        final SettingsModelString m_RESOURCES_FILE_OMNI = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_OMNI,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_omni2.5.hg19.vcf");
        final SettingsModelString m_RESOURCES_FILE_1000G = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1000G,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/1000G_phase1.snps.high_confidence.hg19.vcf");
        final SettingsModelString m_RESOURCES_FILE_DBSNP = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_DBSNP,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/dbsnp_138.hg19.vcf");
        final SettingsModelString m_RESOURCES_FILE_MILLS = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_MILLS,"");// "/storageNGS/ngs1/genomes/mammalian/H_sapiens/hg19_GRCh37/annotation/Mills_and_1000G_gold_standard.indels.hg19.vcf");

        final SettingsModelBoolean m_RESOURCES_BOOLEAN_HAPMAP = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_HAPMAP, true);
        final SettingsModelBoolean m_RESOURCES_BOOLEAN_OMNI = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_OMNI, true);
        final SettingsModelBoolean m_RESOURCES_BOOLEAN_1000G = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1000G, true);
        final SettingsModelBoolean m_RESOURCES_BOOLEAN_DBSNP = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_DBSNP, true);
        final SettingsModelBoolean m_RESOURCES_BOOLEAN_MILLS = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_MILLS, false);

        final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,VQSRNodeModel.DEFAULT_TS_FILTER_SNP,1,100);
    	
        final SettingsModelOptionalString m_OPT_VAR_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_VAR_RECAL,"",false);
        final SettingsModelOptionalString m_OPT_APPLY_RECAL = new SettingsModelOptionalString(VQSRNodeModel.CFGKEY_OPT_APPLY_RECAL,"",false);
    	
    	addPrefPageSetting(m_GATK, IBISKNIMENodesPlugin.GATK);
    	addPrefPageSetting(m_REF_GENOME, IBISKNIMENodesPlugin.REF_GENOME);
    	addPrefPageSetting(m_RESOURCES_FILE_HAPMAP, IBISKNIMENodesPlugin.RES_HAPMAP);
    	addPrefPageSetting(m_RESOURCES_FILE_OMNI, IBISKNIMENodesPlugin.RES_OMNI);
    	addPrefPageSetting(m_RESOURCES_FILE_1000G, IBISKNIMENodesPlugin.RES_1000G_SNPS);
    	addPrefPageSetting(m_RESOURCES_FILE_DBSNP, IBISKNIMENodesPlugin.RES_DBSNP);
    	addPrefPageSetting(m_RESOURCES_FILE_MILLS, IBISKNIMENodesPlugin.RES_MILLS);
    	
    	createNewGroup("General Options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentStringSelection(m_MODE, "Recalibration Mode", "SNP","INDEL"));
    	addDialogComponent(new DialogComponentNumber(m_XMX, "Java Memory per thread in GB",1,4));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("VariantRecalibration");
    	addDialogComponent(new DialogComponentString(m_TRANCHE, "Tranche Levels"));
    	addDialogComponent(new DialogComponentString(m_AN, "Annotation"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_GAUSSIANS, "Gaussians",1,4));
    	addDialogComponent(new DialogComponentNumber(m_NT, "Threads", 1,4));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_VAR_RECAL,"Optional flags"));
    	
    	createNewGroup("Resources");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_HAPMAP, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_HAPMAP, ""));
//    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_HAPMAP, "gatk_vqsr_resource_file1", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_OMNI, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_OMNI, ""));
//    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_OMNI, "gatk_vqsr_resource_file2", 0, ".vcf"));  
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_1000G, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_1000G, ""));
//    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_1000G, "gatk_vqsr_resource_file3", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_DBSNP, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_DBSNP, ""));
//    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_DBSNP, "gatk_vqsr_resource_file4", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(m_RESOURCES_BOOLEAN_MILLS, ""));
    	addDialogComponent(new DialogComponentString(m_RESOURCES_STRING_MILLS, ""));
//    	addDialogComponent(new DialogComponentFileChooser(m_RESOURCES_FILE_MILLS, "gatk_vqsr_resource_file5", 0, ".vcf"));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("ApplyRecalibration");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentNumber(m_TS_FILTER, "TS Filter Level", 0.1,6));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_APPLY_RECAL,"Optional flags"));
    	
//    	createNewTab("VariantRecalibrator");
    	
    	m_RESOURCES_BOOLEAN_HAPMAP.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
//				m_RESOURCES_FILE_HAPMAP.setEnabled(m_RESOURCES_BOOLEAN_HAPMAP.getBooleanValue());
				m_RESOURCES_STRING_HAPMAP.setEnabled(m_RESOURCES_BOOLEAN_HAPMAP.getBooleanValue());
				
			}
		});
    	m_RESOURCES_BOOLEAN_OMNI.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
//				m_RESOURCES_FILE_OMNI.setEnabled(m_RESOURCES_BOOLEAN_OMNI.getBooleanValue());
				m_RESOURCES_STRING_OMNI.setEnabled(m_RESOURCES_BOOLEAN_OMNI.getBooleanValue());	
			}
		});
    	m_RESOURCES_BOOLEAN_1000G.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
//				m_RESOURCES_FILE_1000G.setEnabled(m_RESOURCES_BOOLEAN_1000G.getBooleanValue());
				m_RESOURCES_STRING_1000G.setEnabled(m_RESOURCES_BOOLEAN_1000G.getBooleanValue());
			}
		});
    	m_RESOURCES_BOOLEAN_DBSNP.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
//				m_RESOURCES_FILE_DBSNP.setEnabled(m_RESOURCES_BOOLEAN_DBSNP.getBooleanValue());
				m_RESOURCES_STRING_DBSNP.setEnabled(m_RESOURCES_BOOLEAN_DBSNP.getBooleanValue());	
			}
		});
    	m_RESOURCES_BOOLEAN_MILLS.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
//				m_RESOURCES_FILE_MILLS.setEnabled(m_RESOURCES_BOOLEAN_MILLS.getBooleanValue());
				m_RESOURCES_STRING_MILLS.setEnabled(m_RESOURCES_BOOLEAN_MILLS.getBooleanValue());	
			}
		});
    	
    	m_MODE.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				boolean isSNP = false;
				if(m_MODE.getStringValue().equals("SNP")) {
					isSNP = true;
				}
				m_RESOURCES_BOOLEAN_HAPMAP.setBooleanValue(isSNP);
				m_RESOURCES_STRING_HAPMAP.setEnabled(isSNP);
//				m_RESOURCES_FILE_HAPMAP.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_OMNI.setBooleanValue(isSNP);
				m_RESOURCES_STRING_OMNI.setEnabled(isSNP);
//				m_RESOURCES_FILE_OMNI.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_1000G.setBooleanValue(isSNP);
				m_RESOURCES_STRING_1000G.setEnabled(isSNP);
//				m_RESOURCES_FILE_1000G.setEnabled(isSNP);
				m_RESOURCES_BOOLEAN_DBSNP.setBooleanValue(true);
				m_RESOURCES_STRING_DBSNP.setEnabled(true);
//				m_RESOURCES_FILE_DBSNP.setEnabled(true);
				m_RESOURCES_BOOLEAN_MILLS.setBooleanValue(!isSNP);
				m_RESOURCES_STRING_MILLS.setEnabled(!isSNP);
//				m_RESOURCES_FILE_MILLS.setEnabled(!isSNP);
				
				if(isSNP) {
					m_AN.setStringValue(VQSRNodeModel.DEFAULT_SNP_AN);
					m_TS_FILTER.setDoubleValue(VQSRNodeModel.DEFAULT_TS_FILTER_SNP);
				} else {
					m_AN.setStringValue(VQSRNodeModel.DEFAULT_INDEL_AN);
					m_TS_FILTER.setDoubleValue(VQSRNodeModel.DEFAULT_TS_FILTER_INDEL);
				}
			}
    		
    	});
    }

//	@Override
//    protected void updatePrefs() {
//    	if(usePrefPage.getBooleanValue()) {
//    		String gatkPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("GenomeAnalysisTK.jar");
//    		if(gatkPath != null && !gatkPath.equals("")) {
//    			m_GATK.setStringValue(gatkPath);
//    			m_GATK.setEnabled(false);
//    		} else {
//    			m_GATK.setEnabled(true);
//    		}
//    		String refGenome = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
//    		if(refGenome != null && !refGenome.equals("")) {
//    			m_REF_GENOME.setStringValue(refGenome);
//    			m_REF_GENOME.setEnabled(false);
//    		} else {
//    			m_REF_GENOME.setEnabled(true);
//    		}
//    	} else {
//    		m_GATK.setEnabled(true);
//    		m_REF_GENOME.setEnabled(true);
//    	}
//    }
}

