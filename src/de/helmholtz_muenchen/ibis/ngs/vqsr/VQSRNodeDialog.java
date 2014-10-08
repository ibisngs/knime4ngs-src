package de.helmholtz_muenchen.ibis.ngs.vqsr;

import java.io.File;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

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
public class VQSRNodeDialog extends DefaultNodeSettingsPane {

	
	/**
	 * Models
	 */
    
    private final SettingsModelString m_GATK = new SettingsModelString(VQSRNodeModel.CFGKEY_GATK, "");
    private final SettingsModelString m_REF_GENOME = new SettingsModelString(VQSRNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelString m_INFILE = new SettingsModelString(VQSRNodeModel.CFGKEY_INFILE, "");
    private final SettingsModelString m_MODE = new SettingsModelString(VQSRNodeModel.CFGKEY_MODE, "SNP");
    private final SettingsModelString m_TRANCHE = new SettingsModelString(VQSRNodeModel.CFGKEY_TRANCHE, "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0");
    private final SettingsModelString m_AN = new SettingsModelString(VQSRNodeModel.CFGKEY_AN, "-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum");
    private final SettingsModelIntegerBounded m_NT = new SettingsModelIntegerBounded(VQSRNodeModel.CFGKEY_NT,1,1,Integer.MAX_VALUE);

    private final SettingsModelString m_RESOURCES_STRING_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_1, "Type,known=,training=,truth=,prior=");
    private final SettingsModelString m_RESOURCES_STRING_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_2, "Type,known=,training=,truth=,prior=");
    private final SettingsModelString m_RESOURCES_STRING_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_3, "Type,known=,training=,truth=,prior=");
    private final SettingsModelString m_RESOURCES_STRING_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_STRING_4, "Type,known=,training=,truth=,prior=");
    
    private final SettingsModelString m_RESOURCES_FILE_1 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_1, "");
    private final SettingsModelString m_RESOURCES_FILE_2 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_2, "");
    private final SettingsModelString m_RESOURCES_FILE_3 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_3, "");
    private final SettingsModelString m_RESOURCES_FILE_4 = new SettingsModelString(VQSRNodeModel.CFGKEY_RESOURCES_FILE_4, "");
    
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_1 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_1, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_2 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_2, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_3 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_3, true);
    private final SettingsModelBoolean m_RESOURCES_BOOLEAN_4 = new SettingsModelBoolean(VQSRNodeModel.CFGKEY_RESOURCES_BOOLEAN_4, true);

    private final SettingsModelDoubleBounded m_TS_FILTER = new SettingsModelDoubleBounded(VQSRNodeModel.CFGKEY_TS_FILTER,99.9,1,100);
	
    /**
     * New pane for configuring the VQSR node.
     */
    protected VQSRNodeDialog() {

    
    	createNewGroup("GATK Jar");
    	addDialogComponent(new DialogComponentFileChooser(m_GATK, "gatk_vqsr", 0, ".jar"));

    	createNewGroup("Reference Genome");
    	addDialogComponent(new DialogComponentFileChooser(m_REF_GENOME, "gatk_vqsr_ref_genome", 0, ".fa",".fasta",".txt"));
    	
    	createNewGroup("Infile");
    	addDialogComponent(new DialogComponentFileChooser(m_INFILE, "gatk_vqsr_infile", 0, ".vcf"));
    	
    	createNewGroup("General Options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentStringSelection(m_MODE, "Recalibration Mode", "SNP","INDEL"));
    	addDialogComponent(new DialogComponentString(m_TRANCHE, "Tranche Levels"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentString(m_AN, "Annotation"));
    	addDialogComponent(new DialogComponentNumber(m_NT, "Threads", 1,4));
    	addDialogComponent(new DialogComponentNumber(m_TS_FILTER, "TS Filter Level", 0.1,6));
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
    	
    }
}
