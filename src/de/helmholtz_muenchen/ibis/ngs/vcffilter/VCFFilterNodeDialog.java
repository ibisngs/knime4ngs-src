package de.helmholtz_muenchen.ibis.ngs.vcffilter;

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
import org.knime.core.node.defaultnodesettings.DialogComponentButton;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.defaultnodesettings.SettingsModelStringArray;


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
public class VCFFilterNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the LOFFilter node.
     */
	
	private final SettingsModelString vcf_in = new SettingsModelString(VCFFilterNodeModel.CFGKEY_VCFIN,"-");
	private final SettingsModelString annotation = new SettingsModelString(VCFFilterNodeModel.CFGKEY_ANNOTATION,"");
	private final SettingsModelString vep_script = new SettingsModelString(VCFFilterNodeModel.CFGKEY_VEP_SCRIPT,"");
	
	//LOF definition
	private final SettingsModelString so_term = new SettingsModelString(VCFFilterNodeModel.CFGKEY_SO_TERM,"");
	private final SettingsModelStringArray chosen_terms = new SettingsModelStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST, VCFFilterNodeModel.DEFAULT_TERMS);
	
	private final DialogComponentButton ADD_TERM_BUTTON = new DialogComponentButton("Add selected term");
	
	private static final String NO_SELECTION_MADE = VCFFilterNodeModel.DEFAULT_TERMS[0];
	private final DialogComponentStringListSelection DC_TERM_DISPLAY = new DialogComponentStringListSelection(chosen_terms, "Chosen terms:",NO_SELECTION_MADE);
	
	private final DialogComponentButton REMOVE_TERM_BUTTON = new DialogComponentButton("Remove selected term");
	private final DialogComponentButton RESTORE_DEFAULTS_BUTTON = new DialogComponentButton("Restore default");
	
	//filter
	private final SettingsModelString filter = new SettingsModelString(VCFFilterNodeModel.CFGKEY_FILTER,"");
	private final DialogComponentString DC_FILTER = new DialogComponentString(filter,"Conditions");
	
	HashSet<String> terms = new HashSet<>();
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFFilterNodeModel.class);
	
    protected VCFFilterNodeDialog() {

    	createNewGroup("Select input file");
    	addDialogComponent(new DialogComponentFileChooser(vcf_in, "his_id_vcfin", 0, ".vcf|.vcf.gz"));
    	addDialogComponent(new DialogComponentStringSelection(annotation,"Annotations from",VCFFilterNodeModel.ANNOTATIONS_AVAILABLE));
    	
    	createNewGroup("Path to filter_vep.pl (only VEP)");
    	addDialogComponent(new DialogComponentFileChooser(vep_script, "his_id_vepscript",0, ".pl"));
    	
    	createNewGroup("Select variants");
		addDialogComponent(new DialogComponentStringSelection(so_term,
				"SO terms", VCFFilterNodeModel.SO_TERMS));
		addDialogComponent(ADD_TERM_BUTTON);
		
		DC_TERM_DISPLAY.setVisibleRowCount(5);
		for(String t: VCFFilterNodeModel.DEFAULT_TERMS) {
			terms.add(t);
		}
		DC_TERM_DISPLAY.replaceListItems(terms, NO_SELECTION_MADE);
		addDialogComponent(DC_TERM_DISPLAY);
		
		addDialogComponent(REMOVE_TERM_BUTTON);
		
		RESTORE_DEFAULTS_BUTTON.setToolTipText("Default SO terms comprise variants that are considered as LoF variants.");
		addDialogComponent(RESTORE_DEFAULTS_BUTTON);
		
		createNewGroup("Filter");
		DC_FILTER.setToolTipText("Define additional comma separated filters,e.g.\"LoF is HC,ExAC_AF>0.05");
		DC_FILTER.setSizeComponents(400, 20);
		addDialogComponent(DC_FILTER);

		ADD_TERM_BUTTON.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				addTerm(so_term.getStringValue());
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
				restoreDefault();
			}
		});
    }
    
    public void addTerm(String t) {
    	
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
    
    public void restoreDefault() {
    	terms.clear();
    	for(String t: VCFFilterNodeModel.DEFAULT_TERMS) {
			terms.add(t);
		}
		DC_TERM_DISPLAY.replaceListItems(terms, NO_SELECTION_MADE);
		addDialogComponent(DC_TERM_DISPLAY);
    }
    
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.terms.clear();
    	// check, if data is set
        if (settings.containsKey(VCFFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST))
					this.addTerm(s);
			} catch (InvalidSettingsException e) {
				LOGGER.error(e.getStackTrace());
			}
        }
    }
    
    public void saveAdditionalSettingsTo(NodeSettingsWO settings) {
    	// save the hash to the key FastaSelectorNodeModel.CFGKEY_FILE_LIST
    	settings.addStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST, terms.toArray(new String[terms.size()]));
    }
}

