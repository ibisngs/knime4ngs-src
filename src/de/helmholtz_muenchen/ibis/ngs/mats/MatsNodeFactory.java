package de.helmholtz_muenchen.ibis.ngs.mats;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Mats" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsNodeFactory 
        extends NodeFactory<MatsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MatsNodeModel createNodeModel() {
        return new MatsNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<MatsNodeModel> createNodeView(final int viewIndex,
            final MatsNodeModel nodeModel) {
        return new MatsNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new MatsNodeDialog();
    }
}

