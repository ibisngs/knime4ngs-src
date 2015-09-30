package de.helmholtz_muenchen.ibis.ngs.vepfilter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "LOFFilter" Node.
 * 
 *
 * @author tim.jeske
 */
public class VEPFilterNodeFactory 
        extends NodeFactory<VEPFilterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VEPFilterNodeModel createNodeModel() {
        return new VEPFilterNodeModel();
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
    public NodeView<VEPFilterNodeModel> createNodeView(final int viewIndex,
            final VEPFilterNodeModel nodeModel) {
        return new VEPFilterNodeView(nodeModel);
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
        return new VEPFilterNodeDialog();
    }

}

