package de.helmholtz_muenchen.ibis.ngs.limma;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "Limma" Node.
 * 
 *
 * @author Michael Kluge
 */
public class LimmaNodeView extends NodeView<LimmaNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link LimmaNodeModel})
     */
    protected LimmaNodeView(final LimmaNodeModel nodeModel) {
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

