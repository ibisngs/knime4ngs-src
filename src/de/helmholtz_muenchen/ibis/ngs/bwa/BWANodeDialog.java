package de.helmholtz_muenchen.ibis.ngs.bwa;

import java.io.File;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "BWA" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class BWANodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the BWA node.
     */
    protected BWANodeDialog() {
    	
    	final SettingsModelString readType = new SettingsModelString(BWANodeModel.CFGKEY_READTYPE,"auto-detect");
    	final SettingsModelString refseq = new SettingsModelString(BWANodeModel.CFGKEY_REFSEQFILE,null);
    	final SettingsModelBoolean indexrefseq = new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKINDEX, true);
    	final SettingsModelString readGroup = new SettingsModelString(BWANodeModel.CFGKEY_READGROUP, "");
    	final SettingsModelBoolean readGroupBoolean = new SettingsModelBoolean(BWANodeModel.CFGKEY_READGROUPBOOLEAN, false);
    	readGroup.setEnabled(false);
    	
    	createNewGroup("BWA");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(BWANodeModel.CFGKEY_BWAFILE,"bwa"), "his_id_BWA", 0, ""));
    	createNewGroup("Reference sequence: FastA file (e.g. genome)");
    	addDialogComponent(new DialogComponentFileChooser(refseq, "his1_id_BWA", 0, ""));
    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentBoolean(indexrefseq, "Index reference sequence (Has to be done if index does not exist yet)."));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_BWTINDEX,"BWT-SW"),"Algorithm for constructing BWT index:","BWT-SW","IS"));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_ALNALGO,"BWA-backtrack"),"Algorithm for mapping:","BWA-backtrack","BWA-SW","BWA-MEM"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKCOLORSPACED, false), "Build color-space index."));
    	addDialogComponent(new DialogComponentStringSelection(readType,"Type of reads/ mapping:","auto-detect","single-end","paired-end"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(readGroupBoolean,"Specify Read Group Header:"));
    	addDialogComponent(new DialogComponentString(readGroup,""));
    	setHorizontalPlacement(false);
    	readType.setEnabled(false);
    	refseq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
    	    	String path2refSeq = refseq.getStringValue();
    			File file0 = new File(path2refSeq+".amb");
    			File file1 = new File(path2refSeq+".ann");
    			File file2 = new File(path2refSeq+".bwt");
    			File file3 = new File(path2refSeq+".pac");
    			File file4 = new File(path2refSeq+".sa");
    			if(file0.exists() && file1.exists() && file2.exists() && file3.exists() && file4.exists()){
    				indexrefseq.setBooleanValue(false);
    				indexrefseq.setEnabled(true);
    			} else {
    				indexrefseq.setBooleanValue(true);
    				indexrefseq.setEnabled(false);
    			}
			}
		});
    	
    	readGroupBoolean.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(readGroupBoolean.getBooleanValue()){
						readGroup.setEnabled(true);
					}else{
						readGroup.setEnabled(false);
					}
			}
		});

    	
    }
    
			
}

