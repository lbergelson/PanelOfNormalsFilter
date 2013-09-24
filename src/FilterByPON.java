import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.filters.ClusteredSnps;
import org.broadinstitute.sting.gatk.walkers.filters.FiltrationContext;
import org.broadinstitute.sting.gatk.walkers.filters.FiltrationContextWindow;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.*;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: louisb
 * Date: 9/24/13
 * Time: 1:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class FilterByPON extends RodWalker<Integer, Integer>{

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Input(fullName="panel_of_normals", shortName="panel", doc="Panel of normals", required=false)
    public RodBinding<VariantContext> panel;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;




    @Override
    public void initialize(){
        initializeVcfHeader();
    }

    private void initializeVcfHeader() {
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(this.getToolkit(), rodNames);

        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);


        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        vcfWriter.writeHeader(new VCFHeader(headerLines, vcfSamples));
    }

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Collection<VariantContext> matching = tracker.getValues(panel,ref.getLocus());
        return matching.size();
    }

    @Override
    public Integer reduceInit() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
