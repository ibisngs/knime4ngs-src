package de.helmholtz_muenchen.ibis.ngs.snpeff;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SnpEff" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author
 */
public class SnpEffNodeDialog extends HTExecutorNodeDialog {
    
    public void addToolDialogComponents() {
    	

    	final SettingsModelString snpeff_binary = new SettingsModelString(SnpEffNodeModel.CFGKEY_SNPEFF_BIN, null);
    	final SettingsModelString database = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATABASE, "");
    	final SettingsModelBoolean usebedfile = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_USEBEDFILE, false);
    	final SettingsModelString bed_file = new SettingsModelString(SnpEffNodeModel.CFGKEY_BED_FILE, null);
    	final SettingsModelIntegerBounded memory = new SettingsModelIntegerBounded(SnpEffNodeModel.CFGKEY_MEM, 4, 1, Integer.MAX_VALUE);
        final SettingsModelOptionalString opt_flags = new SettingsModelOptionalString(SnpEffNodeModel.CFGKEY_OPT_FLAGS,"",false);

    	//results filter
    	final SettingsModelBoolean no_downstream = new SettingsModelBoolean(
    			SnpEffNodeModel.CFGKEY_NO_DOWNSTREAM, false);
    	final SettingsModelBoolean no_intergenic = new SettingsModelBoolean(
    			SnpEffNodeModel.CFGKEY_NO_INTERGENIC, false);
    	final SettingsModelBoolean no_intronic = new SettingsModelBoolean(
    			SnpEffNodeModel.CFGKEY_NO_INTRONIC, false);
    	final SettingsModelBoolean no_upstream = new SettingsModelBoolean(
    			SnpEffNodeModel.CFGKEY_NO_UPSTREAM, false);
    	final SettingsModelBoolean no_utr = new SettingsModelBoolean(
    			SnpEffNodeModel.CFGKEY_NO_UTR, false);

    	this.addPrefPageSetting(snpeff_binary, IBISKNIMENodesPlugin.SNPEFF);
    	
    	addDialogComponent(new DialogComponentString(database, "Database name"));
    	addDialogComponent(new DialogComponentBoolean(usebedfile, "Use BED file?"));
    	
    	createNewGroup("Path to BED file");
    	addDialogComponent(new DialogComponentFileChooser(bed_file, "par_3", 0, false, ".bed"));

    	addDialogComponent(new DialogComponentNumber(memory, "Java Memory", 1 ));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentOptionalString(opt_flags,"Optional flags"));
    	
    	createNewTab("Results filter");
    	addDialogComponent(new DialogComponentBoolean(no_downstream, "Don't show DOWNSTREAM changes"));
    	addDialogComponent(new DialogComponentBoolean(no_intergenic, "Don't show INTERGENIC changes"));
    	addDialogComponent(new DialogComponentBoolean(no_intronic, "Don't show INTRONIC changes"));
    	addDialogComponent(new DialogComponentBoolean(no_upstream, "Don't show UPSTREAM changes"));
    	addDialogComponent(new DialogComponentBoolean(no_utr, "Don't show UTR changes"));
    	
    	
    	/*Results filter options change listener*/
    	usebedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usebedfile.getBooleanValue()){
					bed_file.setEnabled(true);
				}
				else{
					bed_file.setEnabled(false);
				}
			}
		});
    }
}

