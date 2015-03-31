package de.helmholtz_muenchen.ibis.ngs.vep;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VEP" Node.
 * 
 *
 * @author tim.jeske
 */
public class VEPNodeFactory 
        extends NodeFactory<VEPNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VEPNodeModel createNodeModel() {
        return new VEPNodeModel();
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
    public NodeView<VEPNodeModel> createNodeView(final int viewIndex,
            final VEPNodeModel nodeModel) {
        return new VEPNodeView(nodeModel);
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
        return new VEPNodeDialog();
    }

}

