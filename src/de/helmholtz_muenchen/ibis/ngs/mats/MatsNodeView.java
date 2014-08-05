package de.helmholtz_muenchen.ibis.ngs.mats;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "Mats" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsNodeView extends NodeView<MatsNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link MatsNodeModel})
     */
    protected MatsNodeView(final MatsNodeModel nodeModel) {
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

