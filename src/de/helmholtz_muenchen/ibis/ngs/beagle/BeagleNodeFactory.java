package de.helmholtz_muenchen.ibis.ngs.beagle;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Beagle" Node.
 * 
 *
 * @author Tanzeem
 */
public class BeagleNodeFactory 
        extends NodeFactory<BeagleNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BeagleNodeModel createNodeModel() {
        return new BeagleNodeModel();
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
    public NodeView<BeagleNodeModel> createNodeView(final int viewIndex,
            final BeagleNodeModel nodeModel) {
        return new BeagleNodeView(nodeModel);
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
        return new BeagleNodeDialog();
    }

}

