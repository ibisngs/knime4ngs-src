package de.helmholtz_muenchen.ibis.ngs.bwa;

import java.io.File;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

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
public class BWANodeDialog extends HTExecutorNodeDialog {

	
	
    /**
     * New pane for configuring the BWA node.
     */
	protected BWANodeDialog() {
	}
	
	public void addToolDialogComponents() {
		
		final SettingsModelString bwa = new SettingsModelString(BWANodeModel.CFGKEY_BWAFILE,"bwa");	
//		final DialogComponentFileChooser dcfc = new DialogComponentFileChooser(bwa,"his_id_BWAPATH", 0);
		final SettingsModelString refseq = new SettingsModelString(BWANodeModel.CFGKEY_REFSEQFILE,"");
		final SettingsModelBoolean indexrefseq = new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKINDEX, true);
		final SettingsModelString readGroup = new SettingsModelString(BWANodeModel.CFGKEY_READGROUP, "@RG\\tID:foo\\tSM:bar");
		final SettingsModelBoolean readGroupBoolean = new SettingsModelBoolean(BWANodeModel.CFGKEY_READGROUPBOOLEAN, false);
		final SettingsModelIntegerBounded ALN_THREADS = new SettingsModelIntegerBounded(BWANodeModel.CFGKEY_THREADS,2, 1, Integer.MAX_VALUE);

    	
    	addPrefPageSetting(bwa, IBISKNIMENodesPlugin.BWA);
    	addPrefPageSetting(refseq, IBISKNIMENodesPlugin.REF_GENOME);
    	readGroup.setEnabled(false);
    	
//    	createNewGroup("BWA");    	
//    	addDialogComponent(dcfc);
    	bwa.setEnabled(false);
//    	createNewGroup("Reference sequence: FastA file (e.g. genome)");
//    	addDialogComponent(new DialogComponentFileChooser(refseq, "his1_id_BWA", 0, ".fa|.fasta"));
//    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentBoolean(indexrefseq, "Index reference sequence (Has to be done if index does not exist yet)."));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_BWTINDEX,"BWT-SW"),"Algorithm for constructing BWT index:","BWT-SW","IS"));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_ALNALGO,"BWA-MEM"),"Algorithm for mapping:","BWA-MEM","BWA-backtrack","BWA-SW"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKCOLORSPACED, false), "Build color-space index."));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(readGroupBoolean,"Specify read group header:"));
    	addDialogComponent(new DialogComponentString(readGroup,""));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentNumber(ALN_THREADS, "Number of threads", 1));
    	
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

