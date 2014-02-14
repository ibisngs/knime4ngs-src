package de.helmholtz_muenchen.ibis.ngs.star;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. * n * nFor details, please see paper: * nA. Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeFactory 
        extends NodeFactory<StarNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public StarNodeModel createNodeModel() {
        return new StarNodeModel();
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
    public NodeView<StarNodeModel> createNodeView(final int viewIndex, final StarNodeModel nodeModel) {
        return new StarNodeView(nodeModel);
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
        return new StarNodeDialog();
    }

}

