package de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;



public abstract class GATKNodeDialog extends HTExecutorNodeDialog{

	
	
    private final SettingsModelString GATK = new SettingsModelString(GATKNodeModel.CFGKEY_GATK_PATH, "");
    private final SettingsModelString REF_GENOME = new SettingsModelString(GATKNodeModel.CFGKEY_REF_GENOME, "");
    private final SettingsModelIntegerBounded m_GATK_MEM = new SettingsModelIntegerBounded(GATKNodeModel.CFGKEY_GATK_MEM, 4, 1, Integer.MAX_VALUE);
    private final SettingsModelOptionalString m_OPT_FLAGS = new SettingsModelOptionalString(GATKNodeModel.CFGKEY_OPT_FLAGS,"",false);

	
    protected GATKNodeDialog() {
    	
    	createNewGroup("Path to GATK jar file");
    	DialogComponentFileChooser gatkf= new DialogComponentFileChooser(GATK, "gatk", JFileChooser.OPEN_DIALOG, false, ".jar");
    	addDialogComponent(gatkf);
    	addDialogComponent(new DialogComponentNumber(m_GATK_MEM, "GATK Memory", 1));
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	addDialogComponent(ref_genome);
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentOptionalString(m_OPT_FLAGS,"Optional flags"));
    	
    	addDialogComponent();
    	
    	usePrefPage.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				GATK.setEnabled(!usePrefPage.getBooleanValue());
				String gatkPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("GenomeAnalysisTK.jar");
		    	if(gatkPath == null) {
		    		gatkPath = "GATK jar not found!";
		    	}
		    	GATK.setStringValue(gatkPath);
			}
    		
    	});

    }
	
    public void onOpen() {
    	super.onOpen();
    	GATK.setEnabled(!usePrefPage.getBooleanValue());
		String gatkPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("GenomeAnalysisTK.jar");
    	if(gatkPath == null) {
    		gatkPath = "GATK jar not found!";
    	}
    	GATK.setStringValue(gatkPath);
    }
	
    protected abstract void addDialogComponent();
    
}
