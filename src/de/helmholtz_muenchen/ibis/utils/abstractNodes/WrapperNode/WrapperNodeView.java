package de.helmholtz_muenchen.ibis.utils.abstractNodes.WrapperNode;


import java.awt.Dimension;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.LogPanel;

/**
 * <code>NodeView</code> which shows STDOUT and STDERR in the view
 * Mostly copied from RNodeView and made RNodeView extend this class in order to avoid code duplication.
 * @author Michael Kluge
 * @author Jonas Zierer
 *
 */
public abstract class WrapperNodeView<T extends WrapperNodeModel> extends NodeView<T> {

    // panel which actually shows the logs
	private final LogPanel PANEL_LOG; 
    
    /**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link ScriptNodeModel}
     */
	protected WrapperNodeView(final T nodeModel) {
        super(nodeModel);
        PANEL_LOG = new LogPanel(nodeModel.getSTDOUT(), nodeModel.getSTDERR());
        PANEL_LOG.setPreferredSize(new Dimension(800, 600));
        setComponent(PANEL_LOG);    
        setViewTitleSuffix("log");
    }

    /** {@inheritDoc} */
    @Override
    protected void modelChanged() {
    	PANEL_LOG.updateView();
    }
    
	@Override
	protected void onClose() {
	
	}

	@Override
	protected void onOpen() {

	}
}
