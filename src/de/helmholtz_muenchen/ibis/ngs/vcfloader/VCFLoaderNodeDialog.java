package de.helmholtz_muenchen.ibis.ngs.vcfloader;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButtonGroup;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentLabel;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VCFLoader" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class VCFLoaderNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring VCFLoader node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
	final SettingsModelString vcf1 = new SettingsModelString(VCFLoaderNodeModel.CFGKEY_VCF1, VCFLoaderNodeModel.DEF_VCF1);
	final SettingsModelString type1 = new SettingsModelString(VCFLoaderNodeModel.CFGKEY_TYPE1, VCFLoaderNodeModel.DEF_TYPE1);
	final SettingsModelBoolean secondfile = new SettingsModelBoolean(VCFLoaderNodeModel.CFGKEY_SECONDFILE, VCFLoaderNodeModel.DEF_SECONDFILE);
	final SettingsModelString vcf2 = new SettingsModelString(VCFLoaderNodeModel.CFGKEY_VCF2, VCFLoaderNodeModel.DEF_VCF2);
	final SettingsModelString type2 = new SettingsModelString(VCFLoaderNodeModel.CFGKEY_TYPE2, VCFLoaderNodeModel.DEF_TYPE2);
	
	
	
    protected VCFLoaderNodeDialog() {
        super();
        
        createNewGroup("VCF file");
        addDialogComponent(new DialogComponentFileChooser(vcf1, "vcfload1", JFileChooser.OPEN_DIALOG, ".vcf"));
        addDialogComponent(new DialogComponentButtonGroup(type1, "Variant type", false, VCFLoaderNodeModel.AVAIL_TYPES , VCFLoaderNodeModel.AVAIL_TYPES));
        createNewGroup("Second VCF file");
        addDialogComponent(new DialogComponentBoolean(secondfile, "Load second VCF file"));
        addDialogComponent(new DialogComponentFileChooser(vcf2, "vcfload2", JFileChooser.OPEN_DIALOG, ".vcf"));
        addDialogComponent(new DialogComponentButtonGroup(type2, "Variant type", false, VCFLoaderNodeModel.AVAIL_TYPES , VCFLoaderNodeModel.AVAIL_TYPES));
        addDialogComponent(new DialogComponentLabel("Variant type has to be different from first file."));
        vcf2.setEnabled(false);
        type2.setEnabled(false);
        
        secondfile.addChangeListener( new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				
				vcf2.setEnabled(secondfile.getBooleanValue());
				type2.setEnabled(secondfile.getBooleanValue());
				
			}
		});
    }
}

