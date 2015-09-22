package de.helmholtz_muenchen.ibis.ngs.vcffilter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

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
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringListSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
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
public class VCFFilterNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VCFFilter node.
     */
	
	private final SettingsModelString vcf_in = new SettingsModelString(VCFFilterNodeModel.CFGKEY_VCFIN,"-");
	private final SettingsModelString vep_script = new SettingsModelString(VCFFilterNodeModel.CFGKEY_VEP_SCRIPT,"");
	private final SettingsModelString vcf_tools = new SettingsModelString(VCFFilterNodeModel.CFGKEY_VCFTOOLS,"-");
	private final SettingsModelBoolean fill_an_ac = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILL_AN_AC,false);

	
	//genotype filter
	private final SettingsModelBoolean filter_by_DP = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_BY_DP,false);
	private final SettingsModelInteger DP_threshold = new SettingsModelInteger(VCFFilterNodeModel.CFGKEY_DP_THRESHOLD,8);
	
	private final SettingsModelBoolean filter_by_GQ = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_BY_GQ,false);
	private final SettingsModelIntegerBounded GQ_threshold = new SettingsModelIntegerBounded(VCFFilterNodeModel.CFGKEY_GQ_THRESHOLD,20,0,99);
	
	//variant filter
	private final SettingsModelBoolean filter_by_GQ_mean = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_BY_GQ_MEAN,false);
	private final SettingsModelIntegerBounded GQ_MEAN_threshold = new SettingsModelIntegerBounded(VCFFilterNodeModel.CFGKEY_GQ_MEAN_THRESHOLD,35,0,99);
	
	private final SettingsModelBoolean filter_pass = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_PASS,false);
	
	private final SettingsModelBoolean filter_by_callRate = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_CALL_RATE,false);
	private final SettingsModelDoubleBounded callRate_threshold = new SettingsModelDoubleBounded(VCFFilterNodeModel.CFGKEY_CALL_RATE_THRESHOLD,0.88,0.0,1.0);
	
	//annotation filter
	private final SettingsModelBoolean filter_annotation = new SettingsModelBoolean(VCFFilterNodeModel.CFGKEY_FILTER_ANNOTATIONS,false);
	private final SettingsModelString annotation = new SettingsModelString(VCFFilterNodeModel.CFGKEY_ANNOTATION,"");
	private final SettingsModelString so_term = new SettingsModelString(VCFFilterNodeModel.CFGKEY_SO_TERM,"");
	private final SettingsModelStringArray chosen_terms = new SettingsModelStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST, VCFFilterNodeModel.DEFAULT_TERMS);
	
	private static final String NO_SELECTION_MADE = VCFFilterNodeModel.DEFAULT_TERMS[0];
	private final DialogComponentStringListSelection DC_TERM_DISPLAY = new DialogComponentStringListSelection(chosen_terms, "Chosen terms:",NO_SELECTION_MADE);
	private final DialogComponentButton ADD_TERM_BUTTON = new DialogComponentButton("Add selected term");
	private final DialogComponentButton REMOVE_TERM_BUTTON = new DialogComponentButton("Remove selected term");
	private final DialogComponentButton RESTORE_DEFAULTS_BUTTON = new DialogComponentButton("Restore default");
	
	private final SettingsModelOptionalString filter = new SettingsModelOptionalString(VCFFilterNodeModel.CFGKEY_FILTER,"",false);
	private final DialogComponentOptionalString DC_FILTER = new DialogComponentOptionalString(filter,"Conditions");
	
	HashSet<String> terms = new HashSet<>();
	
	private static final NodeLogger LOGGER = NodeLogger.getLogger(VCFFilterNodeModel.class);
	
    protected VCFFilterNodeDialog() {

    	createNewGroup("Select input file");
    	addDialogComponent(new DialogComponentFileChooser(vcf_in, "his_id_vcfin", 0, ".vcf|.vcf.gz"));
    	
    	createNewGroup("Path to filter_vep.pl");
    	addDialogComponent(new DialogComponentFileChooser(vep_script, "his_id_vepscript",0, ".pl"));
    	
    	createNewGroup("Path to VCFtools binary");
    	addDialogComponent(new DialogComponentFileChooser(vcf_tools, "his_id_vcftools", 0, ""));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(fill_an_ac,"Restore allele counts?"));

    	
    	//genotype filter
    	createNewTab("Genotype Filter");
		addDialogComponent(new DialogComponentBoolean(filter_by_DP, "Filter genotypes by DP?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(DP_threshold, "DP threshold",1));
		setHorizontalPlacement(false);
		
		addDialogComponent(new DialogComponentBoolean(filter_by_GQ, "Filter genotypes by GQ?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(GQ_threshold, "GQ threshold",1));
		setHorizontalPlacement(false);
		
		//variant filter
		createNewTab("Variant Filter");
		addDialogComponent(new DialogComponentBoolean(filter_pass,"Filter by FILTER is PASS?"));
	
		addDialogComponent(new DialogComponentBoolean(filter_by_callRate,"Filter by call rate?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(callRate_threshold,"Call rate threshold",0.01));
		setHorizontalPlacement(false);
		
		addDialogComponent(new DialogComponentBoolean(filter_by_GQ_mean, "Filter genotypes by mean GQ?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(GQ_MEAN_threshold, "Mean GQ threshold",1));
		setHorizontalPlacement(false);
    	
    	//annotation filter
    	createNewTab("Annotation Filter");
    	addDialogComponent(new DialogComponentBoolean(filter_annotation, "Filter annotations?"));
    	addDialogComponent(new DialogComponentStringSelection(annotation,"Annotations from",VCFFilterNodeModel.ANNOTATIONS_AVAILABLE));
    	
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
		
		usePrefPage.addChangeListener(new ChangeListener (){

			@Override
			public void stateChanged(ChangeEvent e) {
				vcf_tools.setEnabled(!usePrefPage.getBooleanValue());
				vep_script.setEnabled(!usePrefPage.getBooleanValue());
				if(!usePrefPage.getBooleanValue()) return;
		    	
		    	String vcftoolsPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("vcftools");
				if(vcftoolsPath == null) {
					vcftoolsPath = "VCFtools executable not found!";
				}
				vcf_tools.setStringValue(vcftoolsPath);
				
				String filterVep = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("filter_vep.pl");
				if(filterVep == null) {
					filterVep = "filter_vep.pl not found!";
				}
				vep_script.setStringValue(filterVep);
			}
			
		});
		
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
    	for(String t: VCFFilterNodeModel.DEFAULT_TERMS) {
			terms.add(t);
		}
		DC_TERM_DISPLAY.replaceListItems(terms, NO_SELECTION_MADE);
    }
    
    public void loadAdditionalSettingsFrom(NodeSettingsRO settings, DataTableSpec[] specs) throws NotConfigurableException {
    	// clean the old data
    	this.terms.clear();
    	// check, if data is set
        if (settings.containsKey(VCFFilterNodeModel.CFGKEY_TERM_LIST)) {
        	try {
        		// add the values
				for(String s : settings.getStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST)) {
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
    	settings.addStringArray(VCFFilterNodeModel.CFGKEY_TERM_LIST, terms.toArray(new String[terms.size()]));
    }
    
    public void onOpen() {
    	super.onOpen();
    	vcf_tools.setEnabled(!usePrefPage.getBooleanValue());
    	vep_script.setEnabled(!usePrefPage.getBooleanValue());
    	if(!usePrefPage.getBooleanValue()) return;
    	
    	String vcftoolsPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("vcftools");
		if(vcftoolsPath == null) {
			vcftoolsPath = "VCFtools executable not found!";
		}
		vcf_tools.setStringValue(vcftoolsPath);
		
		String filterVep = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("filter_vep.pl");
		if(filterVep == null) {
			filterVep = "filter_vep.pl not found!";
		}
		vep_script.setStringValue(filterVep);
    }
}

