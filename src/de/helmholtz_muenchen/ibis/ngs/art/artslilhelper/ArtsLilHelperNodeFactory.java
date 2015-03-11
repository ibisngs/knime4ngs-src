package de.helmholtz_muenchen.ibis.ngs.art.artslilhelper;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ArtsLilHelper" Node.
 * 
 *
 * @author Tanzeem Haque
 */
public class ArtsLilHelperNodeFactory 
        extends NodeFactory<ArtsLilHelperNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ArtsLilHelperNodeModel createNodeModel() {
        return new ArtsLilHelperNodeModel();
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
    public NodeView<ArtsLilHelperNodeModel> createNodeView(final int viewIndex,
            final ArtsLilHelperNodeModel nodeModel) {
        return new ArtsLilHelperNodeView(nodeModel);
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
        return new ArtsLilHelperNodeDialog();
    }

}

