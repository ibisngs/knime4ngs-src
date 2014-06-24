package de.helmholtz_muenchen.ibis.ngs.matsResultIndexer;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "MatsResultIndexer" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultIndexerNodeView extends NodeView<MatsResultIndexerNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link MatsResultIndexerNodeModel})
     */
    protected MatsResultIndexerNodeView(final MatsResultIndexerNodeModel nodeModel) {
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

