package de.helmholtz_muenchen.ibis.ngs.segemehl;

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
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Segemehl" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class SegemehlNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the Segemehl node.
     */
    protected SegemehlNodeDialog() {

    	final SettingsModelString readType = new SettingsModelString(SegemehlNodeModel.CFGKEY_READTYPE,"single-end");
    	final SettingsModelBoolean clip5adapter = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIP5ADAPTER, false);
    	final SettingsModelBoolean clip3adapter = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIP3ADAPTER, false);
    	final SettingsModelBoolean autoadapter3seq = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_AUTOADAPTER3SEQ, false);
    	final SettingsModelString adapter3seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER3SEQ, "");
    	final SettingsModelString adapter5seq = new SettingsModelString(SegemehlNodeModel.CFGKEY_ADAPTER5SEQ, "");
    	final SettingsModelBoolean clippolya = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CLIPPOLYA, false);
    	final SettingsModelIntegerBounded clippingaccuracy = new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_CLIPPINGACCURACY, 70, 0, 100);
    	final SettingsModelString softhardclipping= new SettingsModelString(SegemehlNodeModel.CFGKEY_SOFTHARDCLIPPING,"Soft (Default)");
    	final SettingsModelBoolean checkBisulfiteMapping = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CHECKSBISULFITEMAPPING, false);
    	final SettingsModelString bisulfiteMappingType = new SettingsModelString(SegemehlNodeModel.CFGKEY_BISULFITEMAPPINGTYPE,"methylC-seq/Lister et al.");
    	final SettingsModelString refseq = new SettingsModelString(SegemehlNodeModel.CFGKEY_REFSEQFILE,null);
    	final SettingsModelBoolean indexrefseq = new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CHECKINDEX, true);
    	autoadapter3seq.setEnabled(false);
    	adapter3seq.setEnabled(false);
    	adapter5seq.setEnabled(false);
    	clippingaccuracy.setEnabled(false);
    	softhardclipping.setEnabled(false);
    	bisulfiteMappingType.setEnabled(false);
    	checkBisulfiteMapping.setEnabled(false);
    	readType.setStringValue("single-end");
    	
    	createNewGroup("Segemehl");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(SegemehlNodeModel.CFGKEY_SEGEMEHLFILE,"bwa"), "his_id_Segemehl", 0, ""));
    	createNewGroup("Reference (e.g. genome) sequence: FastA file.");
    	addDialogComponent(new DialogComponentFileChooser(refseq, "his1_id_Segemehl", 0, ""));
    	createNewGroup("General");
    	addDialogComponent(new DialogComponentBoolean(indexrefseq, "Index reference sequence (Has to be done if index does not exist yet)."));
    	addDialogComponent(new DialogComponentStringSelection(readType,"Type of reads/ mapping:","single-end","paired-end"));
    	addDialogComponent(new DialogComponentNumber(new SettingsModelIntegerBounded(SegemehlNodeModel.CFGKEY_THREADS, 4, 1, 250), "Number of threads/ cores to use:", 1));
    	createNewTab("Further Options");
    	createNewGroup("Alignment parameters");
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(SegemehlNodeModel.CFGKEY_CHECKSPLITREADMAPPING, false), "Use (multiple) split read mapping (e.g. for cDNA reads)"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(checkBisulfiteMapping, "Use bisulfite mapping"));
    	addDialogComponent(new DialogComponentStringSelection(bisulfiteMappingType,"with:","methylC-seq/Lister et al.","bs-seq/Cokus et al. protocol","PAR-CLIP with 4SU","PAR-CLIP with 6SG"));
    	setHorizontalPlacement(false);
    	createNewGroup("Adapter and polyA clipping");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(clip5adapter, "Clip 5' adapters"));
    	addDialogComponent(new DialogComponentString(adapter5seq, "Adapter sequence:", true, 10));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(clip3adapter, "Clip 3' adapters"));
    	addDialogComponent(new DialogComponentString(adapter3seq, "Adapter sequence:", true, 10));
    	addDialogComponent(new DialogComponentBoolean(autoadapter3seq, "Automatic 3’ adapter detection"));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(clippolya, "Automatically clip polyA tails"));
    	addDialogComponent(new DialogComponentNumber(clippingaccuracy, "Clipping accuracy [%]:", 1));
    	addDialogComponent(new DialogComponentStringSelection(softhardclipping,"Type of clipping:","Soft (Default)","Hard"));

    	
    	clip5adapter.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter5seq.setEnabled(clip5adapter.getBooleanValue());
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});

    	clip3adapter.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter3seq.setEnabled(clip3adapter.getBooleanValue());
				autoadapter3seq.setEnabled(clip3adapter.getBooleanValue());
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});
    	
    	autoadapter3seq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				adapter3seq.setEnabled(clip3adapter.getBooleanValue() && !autoadapter3seq.getBooleanValue());
			}
		});
    	
    	clippolya.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				clippingaccuracy.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
				softhardclipping.setEnabled(clippolya.getBooleanValue() || clip5adapter.getBooleanValue() || clip3adapter.getBooleanValue());
			}
		});
    	
    	checkBisulfiteMapping.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				bisulfiteMappingType.setEnabled(checkBisulfiteMapping.getBooleanValue());
			}
		});

    	refseq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
    	    	String path2refSeq = refseq.getStringValue();
    			File file = new File(path2refSeq.substring(0,path2refSeq.lastIndexOf(".")+1)+"idx");    	        
    			if(file.exists()){
    				indexrefseq.setBooleanValue(false);
    				indexrefseq.setEnabled(true);
    			} else {
    				indexrefseq.setBooleanValue(true);
    				indexrefseq.setEnabled(false);
    			}
			}
		});
    	
    	
    }
}

