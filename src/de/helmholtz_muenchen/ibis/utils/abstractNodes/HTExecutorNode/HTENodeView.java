package de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode;

import java.awt.Dimension;

import org.knime.core.node.NodeView;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.LogPanel;

public class HTENodeView<T extends HTExecutorNodeModel> extends NodeView<T>{
    // panel which actually shows the logs
	private final LogPanel PANEL_LOG; 
    
    /**
     * Creates a new view.
     * 
     * @param nodeModel the model class: {@link HTExecutorNodeModel}
     */
	public HTENodeView(final T nodeModel) {
        super(nodeModel);
        PANEL_LOG = new LogPanel(nodeModel.getHTEOUT(), nodeModel.getHTEERR());
        PANEL_LOG.setPreferredSize(new Dimension(800, 600));
        setComponent(PANEL_LOG);    
        setViewTitleSuffix("log");
        setShowNODATALabel(false); // always show view, never show no data
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