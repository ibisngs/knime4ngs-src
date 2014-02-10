package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.edgeRanking;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.util.LinkedHashMap;

import javax.swing.BorderFactory;
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

import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.port.PortObjectSpec;


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
public class GraphicalModelsNodeDialog extends NodeDialogPane {

    /** The node logger for this class. */
    @SuppressWarnings("unused")
	private static final NodeLogger LOGGER = NodeLogger.getLogger(GraphicalModelsNodeDialog.class);

	/** SETTING MODELS */
    private final JList<String> m_typesList;
    private final DefaultListModel<String> m_typesListModel;
    private final JPanel m_typesListPanel;
    private final LinkedHashMap<String, RankerDefinitionPanel> m_typesRankerPanels;
	private final SpinnerNumberModel m_stabSel_sampleNum = new SpinnerNumberModel(10, 1, 10000, 1);
	private final JComboBox<String> m_rankerType;
	

    /**
     * Creates a new binner dialog.
     */
    GraphicalModelsNodeDialog() {
        super();
        
        // INITIALIZE PANEL FOR STABILITY SELECTION PARAMETERS
        final JPanel stabSelPanel = new JPanel(new GridLayout(0,2, 5, 10));
        
        // num samples
        final JLabel stabSel_sampleNum_label     = new JLabel("Number of Samples");
        stabSelPanel.add(stabSel_sampleNum_label);
        final JSpinner stabSel_sampleNum_spinner= new JSpinner(m_stabSel_sampleNum);
        stabSelPanel.add(stabSel_sampleNum_spinner);
        
        // ranking type
        JLabel label_rankerType = new JLabel("Ranking Type");
        stabSelPanel.add(label_rankerType);
        m_rankerType = new JComboBox<String>(GraphicalModelsNodeModel.GRAFO_RANKTYPE);
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
	protected void loadSettingsFrom(final NodeSettingsRO settings, final PortObjectSpec[] specs) throws NotConfigurableException {
		/* NORMAL SETTINGS */
		m_stabSel_sampleNum.setValue(settings.getInt(GraphicalModelsNodeModel.CFGKEY_SS_SSAMPLE_N, 100));
		m_rankerType.setSelectedItem(settings.getString(GraphicalModelsNodeModel.CFGKEY_SS_RANKTYPE , GraphicalModelsNodeModel.GRAFO_RANKTYPE[0]));
		
		/* ADVANCED TYPE/RANKER SETTINGS */
		// delete all old content
        m_typesRankerPanels.clear();
        m_typesListModel.removeAllElements();

		// rankers and params
		String[] types       = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_VAR_TYPES  , new String[0]);
		String[] ranker      = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_RANKER     , new String[0]);
		String[] params      = settings.getStringArray(GraphicalModelsNodeModel.CFGKEY_PARAMS     , new String[0]);

    	
        if(ranker.length != types.length || ranker.length != params.length ){
        	throw(new NotConfigurableException("There has to be one ranker ("+ranker.length+")/ranker params (" + params.length + ") entry per variable type (" + types.length + ")!"));
        }
        
        for(int t=0; t<types.length; t++) {
        	String type = types[t];
        	
        	m_typesListModel.addElement(type);
        	RankerDefinitionPanel p = new RankerDefinitionPanel(type, ranker[t], params[t].replace(GraphicalModelsNodeModel.PARAMS_SEP, "\n"));
			m_typesRankerPanels.put(type, p);
        }
        m_typesList.setSelectedIndex(0);
        getPanel().validate();
        getPanel().repaint();
              
	}


	
	/**
	 * @param settings write intervals to
	 * @throws InvalidSettingsException if a bin name is empty
	 * @see NodeDialogPane#saveSettingsTo(NodeSettingsWO)
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) throws InvalidSettingsException {
		/* NORMAL SETTINGS */
		settings.addInt(GraphicalModelsNodeModel.CFGKEY_SS_SSAMPLE_N, m_stabSel_sampleNum.getNumber().intValue());
		settings.addString(GraphicalModelsNodeModel.CFGKEY_SS_RANKTYPE, (String)this.m_rankerType.getSelectedItem());
		
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
			params[t]   = this.m_typesRankerPanels.get(type).getParams().replace("\n", GraphicalModelsNodeModel.PARAMS_SEP);
		}
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_VAR_TYPES , types);
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_RANKER    , ranker);
		settings.addStringArray(GraphicalModelsNodeModel.CFGKEY_PARAMS    , params);
	}
}

