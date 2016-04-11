package de.helmholtz_muenchen.ibis.ngs.bwa;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "BWA" Node.
 * 
 *
 */
public class BWANodeView extends NodeView<BWANodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link BWANodeModel})
     */
    protected BWANodeView(final BWANodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
    }

}

