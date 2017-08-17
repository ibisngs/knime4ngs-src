/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
import de.helmholtz_muenchen.ibis.utils.datatypes.file.VCFCell;

/**
 * <code>NodeDialog</code> for the "SnpEff" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Sebastian Kopetzky
 */
public class SnpEffNodeDialog extends HTExecutorNodeDialog {
	
	public SnpEffNodeDialog(){
		super(VCFCell.TYPE.getPreferredValueClass(), 0);
	}
    
    public void addToolDialogComponents() {
    	
    	final SettingsModelString snpeff_binary = new SettingsModelString(SnpEffNodeModel.CFGKEY_SNPEFF_BIN, null);
    	final SettingsModelString database = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATABASE, "");
    	final SettingsModelBoolean use_config = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_USE_CONFIG, true);
    	final SettingsModelString database_dir = new SettingsModelString(SnpEffNodeModel.CFGKEY_DATA_DIR, "");
    	final SettingsModelBoolean usebedfile = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_USEBEDFILE, false);
    	final SettingsModelString bed_file = new SettingsModelString(SnpEffNodeModel.CFGKEY_BED_FILE, null);
    	final SettingsModelIntegerBounded memory = new SettingsModelIntegerBounded(SnpEffNodeModel.CFGKEY_MEM, 4, 1, Integer.MAX_VALUE);
        final SettingsModelOptionalString opt_flags = new SettingsModelOptionalString(SnpEffNodeModel.CFGKEY_OPT_FLAGS,"",false);

    	//results filter
    	final SettingsModelBoolean no_downstream = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_NO_DOWNSTREAM, false);
    	final SettingsModelBoolean no_intergenic = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_NO_INTERGENIC, false);
    	final SettingsModelBoolean no_intronic = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_NO_INTRONIC, false);
    	final SettingsModelBoolean no_upstream = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_NO_UPSTREAM, false);
    	final SettingsModelBoolean no_utr = new SettingsModelBoolean(SnpEffNodeModel.CFGKEY_NO_UTR, false);

    	addPrefPageSetting(snpeff_binary, IBISKNIMENodesPlugin.SNPEFF);
    	
    	createNewGroup("Database");
    	addDialogComponent(new DialogComponentString(database, "Name"));
    	addDialogComponent(new DialogComponentBoolean(use_config, "Use database directory specified in config file?"));
    	addDialogComponent(new DialogComponentFileChooser(database_dir, "snpeff_db", 0, true));
    	
    	createNewGroup("Path to BED file");
    	addDialogComponent(new DialogComponentBoolean(usebedfile, "Use BED file?"));
    	addDialogComponent(new DialogComponentFileChooser(bed_file, "par_3", 0, false, ".bed"));

    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentNumber(memory, "Java Memory", 1 ));
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
				bed_file.setEnabled(usebedfile.getBooleanValue());
			}
		});
    	
    	use_config.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent arg0) {
				database_dir.setEnabled(!use_config.getBooleanValue());
			}
    	});
    }
}

