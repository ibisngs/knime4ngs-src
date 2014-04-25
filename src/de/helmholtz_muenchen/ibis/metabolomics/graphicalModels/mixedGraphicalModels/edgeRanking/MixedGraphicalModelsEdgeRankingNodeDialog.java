package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.knime.core.data.DataTableSpec;
import org.knime.core.data.IntValue;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.utils.Global;


/**
 * <code>NodeDialog</code> for the "GraphicalModelsNode" Node.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsEdgeRankingNodeDialog extends NodeDialogPane {

    /** The node logger for this class. */
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(MixedGraphicalModelsEdgeRankingNodeDialog.class);
    
	/** SETTING MODELS */
    private final JList<String> m_typesList;
    private final DefaultListModel<String> m_typesListModel;
    private final JPanel m_typesListPanel;
    private final LinkedHashMap<String, RankerDefinitionPanel> m_typesRankerPanels;
    private final SpinnerNumberModel m_stabSel_sampleNum = new SpinnerNumberModel(MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_N, 1, 10000, 1);
    private final JSpinner stabSel_sampleNum_spinner= new JSpinner(m_stabSel_sampleNum);
    private final SpinnerNumberModel m_stabSel_sampleSize = new SpinnerNumberModel(MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_S, 0.0, 1.0, 0.01);
    private final JComboBox<String> m_randomSeedsSSCol;
    private String[] m_availableCols_ssRSeed = new String[]{MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS};
    private final JComboBox<String> m_rankerType;
	

    /**
     * Creates a new binner dialog.
     */
    protected MixedGraphicalModelsEdgeRankingNodeDialog() {
        super();
        
        // INITIALIZE PANEL FOR STABILITY SELECTION PARAMETERS
        final JPanel stabSelPanel = new JPanel(new GridLayout(0,2, 5, 10));
        
        // num samples
        final JLabel stabSel_sampleNum_label     = new JLabel("Number of Samples");
        stabSelPanel.add(stabSel_sampleNum_label);
        stabSelPanel.add(stabSel_sampleNum_spinner);
        
        // sample Size
        final JLabel stabSel_sampleS_label     = new JLabel("(Relative) Size of Samples");
        stabSelPanel.add(stabSel_sampleS_label);
        final JSpinner stabSel_sampleS_spinner= new JSpinner(m_stabSel_sampleSize);
        stabSelPanel.add(stabSel_sampleS_spinner);
        
        // random Seed
        // Random seeds col if optional input is given
        JLabel label_randomSeed = new JLabel("Random Seed");
        stabSelPanel.add(label_randomSeed);
        m_randomSeedsSSCol = new JComboBox<String>(m_availableCols_ssRSeed);
        stabSelPanel.add(m_randomSeedsSSCol);

        // ranking type
        JLabel label_rankerType = new JLabel("Ranking Type");
        stabSelPanel.add(label_rankerType);
        m_rankerType = new JComboBox<String>(MixedGraphicalModelsEdgeRankingNodeModel.GRAFO_RANKTYPE);
        stabSelPanel.add(m_rankerType);
        
        
        super.addTab(" Stability Selection ", stabSelPanel);
        
        // INITILAIZE PANEL for Ranker Selection
        final JPanel rankerPanel = new JPanel(new GridLayout(1, 1));
        m_typesRankerPanels = new LinkedHashMap<String, RankerDefinitionPanel>();
        
        
        // INITIALIZE TYPES LIST
        m_typesListModel = new DefaultListModel<String>();
        m_typesListModel.addElement("<empty>");
        m_typesList = new JList<String>(m_typesListModel);
        m_typesList.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(final ListSelectionEvent e) {
                currenTypeChanged();
                rankerPanel.validate();
                rankerPanel.repaint();
            }
        });
        m_typesList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        final JScrollPane numScroll = new JScrollPane(m_typesList);
        numScroll.setMinimumSize(new Dimension(200, 155));
        numScroll.setBorder(BorderFactory.createTitledBorder(" Select Type "));

        
        // INITIALIZE 
        m_typesListPanel = new JPanel(new GridLayout(1, 1));
        m_typesListPanel.setBorder(BorderFactory.createTitledBorder(" "));
        m_typesListPanel.setMinimumSize(new Dimension(350, 300));
        m_typesListPanel.setPreferredSize(new Dimension(350, 300));
        JSplitPane split = new JSplitPane(
                JSplitPane.HORIZONTAL_SPLIT, numScroll, m_typesListPanel);
        rankerPanel.add(split);
        super.addTab(" Rankers ", rankerPanel);
    }
    
	private void currenTypeChanged(){
        m_typesListPanel.removeAll();
        String selectedType = m_typesList.getSelectedValue();
        if (selectedType == null) {
            m_typesListPanel.setBorder(BorderFactory.createTitledBorder(" Select Type "));
        } else {
            m_typesListPanel.setBorder(null);
            m_typesListPanel.add(createRankerPanel(selectedType));
        }
	}
	
    private RankerDefinitionPanel createRankerPanel(final String name) {
        RankerDefinitionPanel p;
        if (m_typesRankerPanels.containsKey(name)) {
            p = m_typesRankerPanels.get(name);
        } else {
            p = new RankerDefinitionPanel(name, null, null); // TODO why is this needed?
            m_typesRankerPanels.put(name, p);
        }
        p.validate();
        p.repaint();
        return p;
    }
    
    
    
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadSettingsFrom(final NodeSettingsRO settings, final DataTableSpec[] specs) throws NotConfigurableException {
		HashSet<String> typesCurrent = Global.getUniqueVariableTypes(specs[0]);
		
		/* NORMAL SETTINGS */
		m_stabSel_sampleNum.setValue(settings.getInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_N, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_N));
		m_stabSel_sampleSize.setValue(settings.getDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_SS_SAMPLE_S));
		m_rankerType.setSelectedItem(settings.getString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_RANKTYPE , MixedGraphicalModelsEdgeRankingNodeModel.GRAFO_RANKTYPE[0]));

		if(specs[1] != null){
			DataTableSpec optionalInputSpecs = (DataTableSpec)specs[1];
			this.m_availableCols_ssRSeed = Global.getValidCols(optionalInputSpecs, IntValue.class);
		}else{
			m_availableCols_ssRSeed = new String[]{MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS};
		}
		DefaultComboBoxModel<String> model = new DefaultComboBoxModel<String>( m_availableCols_ssRSeed );
		this.m_randomSeedsSSCol.setModel(model);
		
		if(this.m_availableCols_ssRSeed[0].equals(MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS)){
			this.m_randomSeedsSSCol.setEnabled(false);
			this.stabSel_sampleNum_spinner.setEnabled(true);
		}else{
			this.m_randomSeedsSSCol.setEnabled(true);
			this.stabSel_sampleNum_spinner.setEnabled(false);
			this.m_randomSeedsSSCol.setSelectedItem(settings.getString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANDOM_SEED_COL, MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS));	
		}
		
		/* ADVANCED TYPE/RANKER SETTINGS */
		// delete all old content
        m_typesRankerPanels.clear();
        m_typesListModel.removeAllElements();

		// rankers and params
		String[] types       = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_VAR_TYPES  , new String[0]);
		String[] ranker      = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANKER     , new String[0]);
		String[] params      = settings.getStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARAMS     , new String[0]);

    	
        if(ranker.length != types.length || ranker.length != params.length ){
        	throw(new NotConfigurableException("There has to be one ranker ("+ranker.length+")/ranker params (" + params.length + ") entry per variable type (" + types.length + ")!"));
        }
        
        for(String type: typesCurrent) {
        	int t = Arrays.asList(types).indexOf(type);
        	m_typesListModel.addElement(type);
        	RankerDefinitionPanel p;
        	// ass new
        	if(t<0){
        		p = new RankerDefinitionPanel(type, "", "");
        	}else{	
	        	p = new RankerDefinitionPanel(type, ranker[t], params[t].replace(MixedGraphicalModelsEdgeRankingNodeModel.PARAMS_SEP, "\n"));
        	}
        	m_typesRankerPanels.put(type, p);
        }
        m_typesList.setSelectedIndex(0);
        getPanel().validate();
        getPanel().repaint();
              
	}


	
	/**
	 * @param settings
	 * @throws InvalidSettingsException
	 * @see NodeDialogPane#saveSettingsTo(NodeSettingsWO)
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) throws InvalidSettingsException {
		/* NORMAL SETTINGS */
		settings.addInt(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_N, m_stabSel_sampleNum.getNumber().intValue());
		settings.addDouble(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_SAMPLE_S, m_stabSel_sampleSize.getNumber().doubleValue());
		settings.addString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_SS_RANKTYPE, (String)this.m_rankerType.getSelectedItem());
		settings.addString(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANDOM_SEED_COL, (String)this.m_randomSeedsSSCol.getSelectedItem());

		/* ADVANCED TYPE/RANKER SETTINGS */
		int nTypes = this.m_typesRankerPanels.keySet().size();
		
		// init String arrays
		String[] types      = new String[nTypes];
		String[] ranker     = new String[nTypes];
		String[] params     = new String[nTypes];

		for(int t=0; t<m_typesListModel.size(); t++){
			String type = m_typesListModel.get(t);
			types[t]    = type;
			ranker[t]   = this.m_typesRankerPanels.get(type).getRanker();
			params[t]   = this.m_typesRankerPanels.get(type).getParams().replace("\n", MixedGraphicalModelsEdgeRankingNodeModel.PARAMS_SEP);
		}
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_VAR_TYPES , types);
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_RANKER    , ranker);
		settings.addStringArray(MixedGraphicalModelsEdgeRankingNodeModel.CFGKEY_PARAMS    , params);
	}
}

