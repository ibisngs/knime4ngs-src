package de.helmholtz_muenchen.ibis.ngs.snpsift;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.ngs.snpsift.SnpSiftNodeModel.SnpSiftTool;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SnpSift" Node.
 * 
 * @author Maximilian Hastreiter
 */
public class SnpSiftNodeDialog extends HTExecutorNodeDialog {
    
    public void addToolDialogComponents() {

    	final SettingsModelString snpsift_bin = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_SNPSIFT_BIN,"");
    	final SettingsModelString method = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_METHOD,SnpSiftNodeModel.DEF_METHOD);
    	
    	/**Filter**/
    	final SettingsModelString filterstring = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_FILTERSTRING,"");
    	
    	/**Annotate**/
    	final SettingsModelOptionalString anninfo = new SettingsModelOptionalString(
    			SnpSiftNodeModel.CFGKEY_ANNINFO,SnpSiftNodeModel.DEF_FIELDS,false);
    	final SettingsModelBoolean annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
    	final SettingsModelString annvcfdb = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
    	final SettingsModelOptionalString ann_opt = new SettingsModelOptionalString(
    			SnpSiftNodeModel.CFGKEY_ANN_OPT, "", false);
    	
    	/**Intervals**/
    	final SettingsModelString interbed = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_INTERBED,"");
    	final SettingsModelBoolean interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

    	/**dbNSFP**/
    	final SettingsModelString dbnsfp = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
    	final SettingsModelOptionalString dbnsfpfields = new SettingsModelOptionalString(
    			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,SnpSiftNodeModel.DEF_FIELDS,false);
    	final SettingsModelBoolean dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);
    	final SettingsModelOptionalString dbnsfp_opt = new SettingsModelOptionalString(SnpSiftNodeModel.CFGKEY_DBNSFP_OPT, "", false);

    	
    	/**Other**/
    	final SettingsModelString other_cmd = new SettingsModelString(SnpSiftNodeModel.CFGKEY_OTHER_CMD,"");
    	final SettingsModelString other_out = new SettingsModelString(SnpSiftNodeModel.CFGKEY_OTHER_OUT,"");
    	
    	this.addPrefPageSetting(snpsift_bin, IBISKNIMENodesPlugin.SNPSIFT);
   
    	addDialogComponent(new DialogComponentStringSelection(method, "Select tool",SnpSiftNodeModel.NAME2TOOL.keySet()));
    	
    	createNewGroup("Individual Command");
    	addDialogComponent(new DialogComponentString(other_cmd, "Your command"));
    	addDialogComponent(new DialogComponentFileChooser(other_out, "Output file"));
    	
    	createNewTab("Filter");
    	addDialogComponent(new DialogComponentString(filterstring, "Filter criteria (SnpSift Syntax)"));
    	
    	createNewTab("Annotate");
    	createNewGroup("VCF file providing annotations");
    	addDialogComponent(new DialogComponentFileChooser(annvcfdb, "par_3", 0, false, ".vcf"));
    	addDialogComponent(new DialogComponentBoolean(annid, "Do not annotate INFO fields"));
    	addDialogComponent(new DialogComponentOptionalString(anninfo, "List of INFO fields"));
    	addDialogComponent(new DialogComponentOptionalString(ann_opt, "Futher annotate flags"));
    	
    	createNewTab("Intervals");
    	createNewGroup("Specify interval file");
    	addDialogComponent(new DialogComponentFileChooser(interbed, "par_4", 0, false, ".bed",".txt"));
    	addDialogComponent(new DialogComponentBoolean(interx, "Exclude entries in intervals"));    	
    	
    	createNewTab("dbNSFP");
    	createNewGroup("Specify dbNSFP database");
    	addDialogComponent(new DialogComponentFileChooser(dbnsfp, "par_5", 0, false,".txt.gz"));
    	addDialogComponent(new DialogComponentBoolean(dbnsfpfieldsall, "Include empty values"));   
    	addDialogComponent(new DialogComponentOptionalString(dbnsfpfields, "List of field names"));
    	addDialogComponent(new DialogComponentOptionalString(dbnsfp_opt, "Further dbnsfp flags"));
    	
    	
    	method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				
				other_cmd.setEnabled(false);
				other_out.setEnabled(false);
				
				setEnabled(false, "Filter");
				setEnabled(false, "Annotate");
				setEnabled(false, "Intervals");
				setEnabled(false, "dbNSFP");
				
				SnpSiftTool tool = SnpSiftNodeModel.NAME2TOOL.get(method.getStringValue());
				switch(tool) {
				case FILTER:
					setEnabled(true, "Filter");
					break;
				case ANNOTATE:
					setEnabled(true, "Annotate");
					break;
				case INTERVALS:
					setEnabled(true, "Intervals");
					break;
				case DBNSFP:
					setEnabled(true, "dbNSFP");
					break;
				case OTHER:
					other_cmd.setEnabled(true);
					other_out.setEnabled(true);
					break;
				default:
					break;
				}		
			}
		});
    }
}

