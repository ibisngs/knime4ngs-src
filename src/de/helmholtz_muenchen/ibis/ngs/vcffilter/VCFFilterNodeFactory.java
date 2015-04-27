package de.helmholtz_muenchen.ibis.ngs.vcffilter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "LOFFilter" Node.
 * 
 *
 * @author tim.jeske
 */
public class VCFFilterNodeFactory 
        extends NodeFactory<VCFFilterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFFilterNodeModel createNodeModel() {
        return new VCFFilterNodeModel();
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
    public NodeView<VCFFilterNodeModel> createNodeView(final int viewIndex,
            final VCFFilterNodeModel nodeModel) {
        return new VCFFilterNodeView(nodeModel);
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
        return new VCFFilterNodeDialog();
    }

}

