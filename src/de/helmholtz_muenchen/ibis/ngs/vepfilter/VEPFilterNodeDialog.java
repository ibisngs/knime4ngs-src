package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashSet;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "VEPFilter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class VEPFilterNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VCFFilter node.
     */
	
	private final SettingsModelString vep_script = new SettingsModelString(VEPFilterNodeModel.CFGKEY_VEP_SCRIPT,"");
	
	//annotation filter
	private final SettingsModelString so_term = new SettingsModelString(VEPFilterNodeModel.CFGKEY_SO_TERM,"");
	private final SettingsModelStringArray chosen_terms = new SettingsModelStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST, VEPFilterNodeModel.DEFAULT_TERMS);
	
	private static final String NO_SELECTION_MADE = VEPFilterNodeModel.DEFAULT_TERMS[0];
	private final DialogComponentStringListSelection DC_TERM_DISPLAY = new DialogComponentStringListSelection(chosen_terms, "Chosen terms:",NO_SELECTION_MADE);
	private final DialogComponentButton ADD_TERM_BUTTON = new DialogComponentButton("Add selected term");
	private final DialogComponentButton REMOVE_TERM_BUTTON = new DialogComponentButton("Remove selected term");
	private final DialogComponentButton RESTORE_DEFAULTS_BUTTON = new DialogComponentButton("Restore default");
	
	private final SettingsModelOptionalString filter = new SettingsModelOptionalString(VEPFilterNodeModel.CFGKEY_FILTER,"",false);
	private final DialogComponentOptionalString DC_FILTER = new DialogComponentOptionalString(filter,"Conditions");
	
	private final SettingsModelString outfolder = new SettingsModelString(VEPFilterNodeModel.CFGKEY_OUTFOLDER,"");
	private final SettingsModelBoolean overwrite = new SettingsModelBoolean(VEPFilterNodeModel.CFGKEY_OVERWRITE,false);
	
	HashSet<String> terms = new HashSet<>();
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VEPFilterNodeModel.class);
	
    protected VEPFilterNodeDialog() {
    	
    	createNewGroup("Path to filter_vep.pl");
    	addDialogComponent(new DialogComponentFileChooser(vep_script, "his_id_vepscript",0, ".pl"));
    	
    	createNewGroup("Select variants");
		addDialogComponent(new DialogComponentStringSelection(so_term,
				"SO terms", VEPFilterNodeModel.SO_TERMS));
		addDialogComponent(ADD_TERM_BUTTON);
		
		
		
		DC_TERM_DISPLAY.setVisibleRowCount(5);
		for(String t: VEPFilterNodeModel.DEFAULT_TERMS) {
			terms.add(t);
		}
		DC_TERM_DISPLAY.replaceListItems(terms, NO_SELECTION_MADE);
		addDialogComponent(DC_TERM_DISPLAY);
		
		addDialogComponent(REMOVE_TERM_BUTTON);
		
		RESTORE_DEFAULTS_BUTTON.setToolTipText("Default SO terms comprise variants that are considered as LoF variants.");
		addDialogComponent(RESTORE_DEFAULTS_BUTTON);
		
		createNewGroup("Filter");
		DC_FILTER.setToolTipText("Define additional comma separated filters,e.g.\"LoF is HC,ExAC_AF>0.05");
//		DC_FILTER.setSizeComponents(400, 20);
		addDialogComponent(DC_FILTER);

		ADD_TERM_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				addSOTerm(so_term.getStringValue());
			}
		});
		
		REMOVE_TERM_BUTTON.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				removeTerm(chosen_terms.getStringArrayValue());
			}
		});
		
		RESTORE_DEFAULTS_BUTTON.addActionListener(new ActionListener () {
			@Override
			public void actionPerformed(ActionEvent e) {
				restoreTermDefaults();
			}
		});
		
		createNewGroup("Folder for output files");
    	addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_VEP_Filtered_OUT", 0, true));
    	addDialogComponent(new DialogComponentBoolean(overwrite, "Overwrite, if output files exist?"));
    }
    
    public void addSOTerm(String t) {
    	
    	if(terms.contains(t)) return;
    	
    	terms.add(t);
    	DC_TERM_DISPLAY.replaceListItems(terms, new String[0]);
    }
    
    public void removeTerm(String [] terms) {
    	for(String t:terms) {
    		this.terms.remove(t);
    	}
    	
    	if(this.terms.size()==0) {
    		ArrayList<String> empty = new ArrayList<String>();
			empty.add("--no term selected--");
			DC_TERM_DISPLAY.replaceListItems(empty, "--no term selected--");
    	} else {
    		DC_TERM_DISPLAY.replaceListItems(this.terms, (String[]) null);
    	}
    }
    
    public void restoreTermDefaults() {
    	terms.clear();
    	for(String t: VEPFilterNodeModel.DEFAULT_TERMS) {
			terms.add(t);
		}
		DC_TERM_DISPLAY.replaceListItems(terms, NO_SELECTION_MADE);
    }
    
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.terms.clear();
    	// check, if data is set
        if (settings.containsKey(VEPFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST)) {
					this.addSOTerm(s);
				}
				if(this.terms.size()==0) {
					ArrayList<String> empty = new ArrayList<String>();
					empty.add("--no term selected--");
					DC_TERM_DISPLAY.replaceListItems(empty, "--no term selected--");
				}
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }
    
    public void saveAdditionalSettingsTo(NodeSettingsWO settings) {
    	settings.addStringArray(VEPFilterNodeModel.CFGKEY_TERM_LIST, terms.toArray(new String[terms.size()]));
    }
    
	@Override
	protected void updatePrefs() {
		if(usePrefPage.getBooleanValue()) {
			String filterVep = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("filter_vep.pl");
			if(filterVep != null && !filterVep.equals("")) {
				vep_script.setStringValue(filterVep);
				vep_script.setEnabled(false);
			} else {
				vep_script.setEnabled(true);
			}
		} else {
			vep_script.setEnabled(true);
		}
	}
}