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
    			SnpSiftNodeModel.CFGKEY_ANNINFO,"",false);
    	final SettingsModelBoolean annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
    	final SettingsModelString annvcfdb = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
    	final SettingsModelOptionalString ann_opt = new SettingsModelOptionalString(
    			SnpSiftNodeModel.CFGKEY_ANN_OPT, "", false);
    	
    	/**Intervals**/
    	final SettingsModelString interbed = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_INTERBED,"");
    	final SettingsModelBoolean interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

    	/**dbnsfp**/
    	final SettingsModelString dbnsfp = new SettingsModelString(
    			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
    	final SettingsModelOptionalString dbnsfpfields = new SettingsModelOptionalString(
    			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,"",false);
    	final SettingsModelBoolean dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);
    	
    	this.addPrefPageSetting(snpsift_bin, IBISKNIMENodesPlugin.SNPSIFT);
   
    	addDialogComponent(new DialogComponentStringSelection(method, "Select tool",SnpSiftNodeModel.NAME2TOOL.keySet()));
    	addDialogComponent(new DialogComponentString(filterstring, "Filter criteria (SnpSift Syntax)"));
		annvcfdb.setEnabled(false);
    	
    	createNewTab("Annotate");
    	createNewGroup("VCF file providing annotations");
    	addDialogComponent(new DialogComponentFileChooser(annvcfdb, "par_3", 0, false, ".vcf"));
//    	closeCurrentGroup();
    	addDialogComponent(new DialogComponentBoolean(annid, "Do not annotate INFO fields"));
    	addDialogComponent(new DialogComponentOptionalString(anninfo, "Annotate using a (comma-separated) list of info fields"));
    	addDialogComponent(new DialogComponentOptionalString(ann_opt, "Futher annotate flags"));
    	
    	createNewTab("Intervals");
    	createNewGroup("Specify interval file");
    	addDialogComponent(new DialogComponentFileChooser(interbed, "par_4", 0, false, ".bed",".txt"));
    	addDialogComponent(new DialogComponentBoolean(interx, "Exclude VCF entries in intervals"));    	
    	
    	createNewTab("Annotate with dbnsfp");
    	createNewGroup("Specify dbnsfp database");
    	addDialogComponent(new DialogComponentFileChooser(dbnsfp, "par_5", 0, false,".txt"));
    	addDialogComponent(new DialogComponentBoolean(dbnsfpfieldsall, "Annotate all fields, even if the database has an empty value"));   
    	addDialogComponent(new DialogComponentOptionalString(dbnsfpfields, "Specify which fields are used for annotation by a comma separated list of field names"));

    	dbnsfpfieldsall.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(dbnsfpfieldsall.getBooleanValue()){
						dbnsfpfields.setEnabled(false);
						dbnsfpfields.setIsActive(false);
					}else{
						dbnsfpfields.setEnabled(true);
					}
			}
		});
    	
    	method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				
				filterstring.setEnabled(false);
				setEnabled(false, "Annotate");
//				setEnabled(true, "Intervals");
//				setEnabled(false, "Annotate with dbnsfp");
//				setEnabled(true, method.getStringValue());
				
				SnpSiftTool tool = SnpSiftNodeModel.NAME2TOOL.get(method.getStringValue());
				switch(tool) {
				case FILTER:
					filterstring.setEnabled(true);
					break;
				case ANNOTATE:
					setEnabled(true, "Annotate");
				default:
					break;
				}		
			}
		});
    }
}

