package de.helmholtz_muenchen.ibis.ngs.loffilter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "LOFFilter" Node.
 * 
 *
 * @author tim.jeske
 */
public class LOFFilterNodeFactory 
        extends NodeFactory<LOFFilterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public LOFFilterNodeModel createNodeModel() {
        return new LOFFilterNodeModel();
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
    public NodeView<LOFFilterNodeModel> createNodeView(final int viewIndex,
            final LOFFilterNodeModel nodeModel) {
        return new LOFFilterNodeView(nodeModel);
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
        return new LOFFilterNodeDialog();
    }

}

