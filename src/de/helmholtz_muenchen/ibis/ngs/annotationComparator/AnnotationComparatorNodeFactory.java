package de.helmholtz_muenchen.ibis.ngs.annotationComparator;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "AnnotationComparator" Node.
 * 
 *
 * @author tim.jeske
 */
public class AnnotationComparatorNodeFactory 
        extends NodeFactory<AnnotationComparatorNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public AnnotationComparatorNodeModel createNodeModel() {
        return new AnnotationComparatorNodeModel();
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
    public NodeView<AnnotationComparatorNodeModel> createNodeView(final int viewIndex,
            final AnnotationComparatorNodeModel nodeModel) {
        return new AnnotationComparatorNodeView(nodeModel);
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
        return new AnnotationComparatorNodeDialog();
    }

}

