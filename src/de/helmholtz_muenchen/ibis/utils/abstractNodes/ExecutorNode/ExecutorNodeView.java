package de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode;


import java.awt.Dimension;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ScriptNode.LogPanel;

/**
 * <code>ExecutorNodeView</code> which shows STDOUT and STDERR in the view
 * 
 * @author Jonas Zierer
 * @author Michael Kluge
 *
 */
public abstract class ExecutorNodeView<T extends ExecutorNodeModel> extends NodeView<T> {

    // panel which actually shows the logs
	private final LogPanel PANEL_LOG; 
    
    /**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link ScriptNodeModel}
     */
	protected ExecutorNodeView(final T nodeModel) {
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
