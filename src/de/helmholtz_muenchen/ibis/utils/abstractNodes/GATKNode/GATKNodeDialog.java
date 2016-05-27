package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.depthofcoverage.DepthOfCoverageNodeModel;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;



public abstract class GATKNodeDialog extends HTExecutorNodeDialog{

//    private final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
//    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
//    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
//	  private final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
//	  private final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
//    private final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_OPT_FLAGS,"",false);

    protected GATKNodeDialog() {
    }
    
    public void addToolDialogComponents() {
    	addDialogComponent();
    	
    	if(optionsTabIsEmpty()) {
			setDefaultTabTitle("GATK");
		} else {
			createNewTab("GATK");
		}
    	
    	final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
        final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
        final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
    	final SettingsModelString m_path2bed = new SettingsModelOptionalString(DepthOfCoverageNodeModel.CFGKEY_PATH2BED,"",true);
    	final SettingsModelBoolean m_bed_file_check = new SettingsModelBoolean(DepthOfCoverageNodeModel.CFGKEY_BED_FILE_CHECKBOX,false);
        final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_OPT_FLAGS,"",false);
    	
    	addPrefPageSetting(REF_GENOME, IBISKNIMENodesPlugin.REF_GENOME);
    	addPrefPageSetting(GATK, IBISKNIMENodesPlugin.GATK);
    	
    	addDialogComponent(new DialogComponentNumber(m_GATK_MEM, "GATK Memory", 1));
    	addDialogComponent(new DialogComponentBoolean(m_bed_file_check,"Use BED file?"));
    	m_path2bed.setEnabled(false);
    	
    	createNewGroup("Path to BED file");
    	addDialogComponent(new DialogComponentFileChooser(m_path2bed, "his_id_GATK_DoC", 0, ".bed"));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_FLAGS,"Optional flags"));
    	
    	
    	m_bed_file_check.addChangeListener(new ChangeListener () {
    		
    		@Override
			public void stateChanged(ChangeEvent e) {
				m_path2bed.setEnabled(m_bed_file_check.getBooleanValue());
			}
    	});

    }
	
    protected abstract void addDialogComponent();
    
}
