package de.helmholtz_muenchen.ibis.plotting.scatterplot;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.data.DoubleValue;
import org.knime.core.data.IntValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RPlottingNodeDialog;

/**
 * <code>NodeDialog</code> for the "Boxplot" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class ScatterplotNodeDialog extends RPlottingNodeDialog {

    // models
    protected final SettingsModelBoolean m_matrix = new SettingsModelBoolean(ScatterplotNodeModel.CFGKEY_MATRIX, false);
    private final SettingsModelFilterString m_matrix_columns   = new SettingsModelFilterString(ScatterplotNodeModel.CFGKEY_MATRIX_COLUMNS);
    protected final SettingsModelDoubleBounded m_alpha = new SettingsModelDoubleBounded(ScatterplotNodeModel.CFGKEY_ALPHA, 0.8, 0.0, 1.0);
    protected final SettingsModelIntegerBounded m_pointsize = new SettingsModelIntegerBounded(ScatterplotNodeModel.CFGKEY_POINTSIZE, 1, 1, 20);

	
    /**
     * New pane for configuring the Boxplot node.
     */
	@SuppressWarnings("unchecked")
	protected ScatterplotNodeDialog() {
        super();        
        
        this.createNewTabAt("Scatterplot Options", 0);
        addDialogComponent(new DialogComponentNumber(
        		m_alpha,
        		"Alpha", 
        		0.1, 3));
        
        addDialogComponent(new DialogComponentNumber(
        		m_pointsize,
        		"Pointsize", 
        		1, 3));
        
        this.createNewGroup("Scatterplot Matrix");
        addDialogComponent(new DialogComponentBoolean(
        		m_matrix, 
        		"Scatterplot Matrix"));
        m_matrix.addChangeListener(new ChangeListener() {
      	    public void stateChanged(final ChangeEvent e) {
      	    	if(m_matrix.getBooleanValue()){
      	    		m_facet.setEnabled(false);
      	    		m_matrix_columns.setEnabled(true);
      	    		m_col_x.setEnabled(false);
      	    		m_col_y.setEnabled(false);
      	    	}else{
      	    		m_facet.setEnabled(true);
      	    		m_matrix_columns.setEnabled(false);
      	    		m_col_x.setEnabled(true);
      	    		m_col_y.setEnabled(true);
      	    	}
      	    }
      	});
		// Variables
		addDialogComponent(new DialogComponentColumnFilter(
				m_matrix_columns,
				0, true, IntValue.class, DoubleValue.class));
		this.closeCurrentGroup();
  
    }
}

