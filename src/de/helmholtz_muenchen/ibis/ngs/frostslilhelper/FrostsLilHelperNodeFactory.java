package de.helmholtz_muenchen.ibis.ngs.frostslilhelper;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FrostsLilHelper" Node.
 * 
 *
 * @author 
 */
public class FrostsLilHelperNodeFactory 
        extends NodeFactory<FrostsLilHelperNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FrostsLilHelperNodeModel createNodeModel() {
        return new FrostsLilHelperNodeModel();
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
    public NodeView<FrostsLilHelperNodeModel> createNodeView(final int viewIndex,
            final FrostsLilHelperNodeModel nodeModel) {
        return new FrostsLilHelperNodeView(nodeModel);
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
        return new FrostsLilHelperNodeDialog();
    }

}

