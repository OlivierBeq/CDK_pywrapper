package CDKdesc;

import org.apache.commons.cli.*;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.fingerprint.*;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.AtomTypeAwareSaturationChecker;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Main {
    public static void main(String[] args) throws Exception {
        // Set switches to be used
        CommandLineParser parser = new DefaultParser();
        Options options = new Options();

        options.addOption("f", "fingerprint", true, "Calculate a specific fingerprint. " +
                "\tMust be one of {FP, ExtFP, EStateFP, GraphFP, MACCSFP, PubchemFP,\n" +
                "\tSubFP, KRFP, AP2DFP, HybridFP, LingoFP, SPFP, SigFP, CircFP}");
        options.addOption("nBits", true, "Number of bits of FP and GraphFP fingerprints (default: 1024)");
        options.addOption("depth", true, "Search depth of FP and GraphFP fingerprints (default: 6)");
        options.addOption("i", "input", true, "Input v2000 SD file");
        options.addOption("h", "help", false, "Shows this Help");

        try {
            // Parse switches
            CommandLine commandLine = parser.parse(options, args);
            // Obtain a descriptor engine
            DescriptorEngine descriptorEngine = new DescriptorEngine(IMolecularDescriptor.class,
                    SilentChemObjectBuilder.getInstance());
            List<IDescriptor> descriptors = descriptorEngine.getDescriptorInstances();
            descriptors.add(new JPlogPDescriptor());
            // Fingerprint names, to be chosen from
            List<String> fp_names = new ArrayList<>(Arrays.asList("FP", "ExtFP", "EStateFP", "GraphFP", "MACCSFP",
                    "PubchemFP", "SubFP", "KRFP", "AP2DFP", "HybridFP", "LingoFP", "SPFP", "SigFP", "CircFP"));
            // Display help
            if (commandLine.hasOption("help")) {
                new HelpFormatter().printHelp("java -jar CDKdesc.jar", options);
            } else if (!commandLine.hasOption("input")){
                // No input given
                throw new Exception("Input V2000 SD file must be provided.");
            } else if (commandLine.hasOption("fingerprint")) {
                // Obtain type of fingerprint
                String fp_value = commandLine.getOptionValue("fingerprint");
                if (!fp_names.contains(fp_value)){
                    // Not supported
                    throw new Exception("Fingerprint type " + fp_value + " is not available.");
                }
                int size = 1024;
                // Get size of fingerprint
                if (commandLine.hasOption("nBits")) {
                    size = Integer.parseInt(commandLine.getOptionValue("nBits"));
                }
                // Same with search depth
                int searchDepth = 6;
                if (commandLine.hasOption("searchDepth")) {
                    searchDepth = Integer.parseInt(commandLine.getOptionValue("searchDepth"));
                }
                // Instantiate fingerprint calculator
                IFingerprinter fp = switch (fp_value) {
                    case "FP" -> new Fingerprinter(size, searchDepth);
                    case "ExtFP" -> new ExtendedFingerprinter(size, searchDepth);
                    case "EStateFP" -> new EStateFingerprinter();
                    case "GraphFP" -> new GraphOnlyFingerprinter(size, searchDepth);
                    case "MACCSFP" -> new MACCSFingerprinter();
                    case "PubchemFP" -> new PubchemFingerprinter(SilentChemObjectBuilder.getInstance());
                    case "SubFP" -> new SubstructureFingerprinter();
                    case "KRFP" -> new KlekotaRothFingerprinter();
                    case "AP2DFP" -> new AtomPairs2DFingerprinter();
                    case "HybridFP" -> new HybridizationFingerprinter(size, searchDepth);
                    case "LingoFP" -> new LingoFingerprinter(searchDepth);
                    case "SPFP" -> new ShortestPathFingerprinter(size);
                    case "SigFP" -> new SignatureFingerprinter(searchDepth);
                    case "CircFP" -> new CircularFingerprinter();
                    default ->
                        // Default to eCDKFingerprinter
                            throw new Exception("Fingerprint type " + fp_value + " is not implemented.");
                };
                // Obtain FP names and print
                List<String> value_names = new ArrayList<>(); //
                boolean obtained_names = false;

                // Calculate fingerprint values
                // Read content of file
                try (FileInputStream fis = new FileInputStream(commandLine.getOptionValue("input"))) {
                    IteratingSDFReader supplier = new IteratingSDFReader(fis, DefaultChemObjectBuilder.getInstance());
                    // Iterate over molecules
                    int idx_mol = 0;
                    while (supplier.hasNext()) {
                        IAtomContainer molecule = supplier.next();
                        // Do perception (atom types & aromaticity)
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
                                Cycles.vertexShort());
                        aromaticity.apply(molecule);
                        // Calculate fingerprint
                        IBitFingerprint fp_values = fp.getBitFingerprint(molecule);
                        int fp_size;
                        if (fp_values instanceof BitSetFingerprint){
                            BitSet fp_bitset = fp_values.asBitSet();
                            // Obtain FP size
                            fp_size = fp.getSize();
                            if (fp_size < 0){
                                fp_size = (int)fp_values.size();
                            }
                            String[] bits = new String[fp_size];
                            Arrays.fill(bits, "0");
                            // Obtain first set bit
                            int idx = fp_bitset.nextSetBit(0);
                            while (idx >= 0) {
                                bits[idx] = "1";
                                idx = fp_bitset.nextSetBit(idx + 1);
                            }
                            if (!obtained_names) {
                                value_names = (IntStream.rangeClosed(1, fp_size)
                                        .mapToObj(x -> fp_value + "_" + x)
                                        .toList());
                                obtained_names = true;
                                System.out.println(String.join(" ", value_names));
                            }
                            // Print values to sdtout
                            System.out.println(String.join(" ", bits));
                        } else { // Signature is a Hashmap
                            Map<String, Integer> raw_fp = fp.getRawFingerprint(molecule);
                            String[] bits = new String[raw_fp.size()];
                            int index = 0;
                            for (Map.Entry<String, Integer> entry : raw_fp.entrySet()){
                                bits[index] = "\"" + entry.getKey() + "\": " + entry.getValue().toString();
                                index += 1;
                            }
                            System.out.println(String.valueOf(idx_mol) + ": {" + String.join(", ", bits) + "},");
                            idx_mol +=1;
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            } else {
                // Calculate descriptor values
                List<String> value_names = new ArrayList<>();
                boolean obtained_names = false;
                // Read content of file
                try (FileInputStream fis = new FileInputStream(commandLine.getOptionValue("input"))) {
                    IteratingSDFReader supplier = new IteratingSDFReader(fis, DefaultChemObjectBuilder.getInstance());
                    // Iterate over molecules
                    while (supplier.hasNext()) {
                        IAtomContainer molecule = supplier.next();
                        // Do perception (atom types & aromaticity)
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                        Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(),
                                Cycles.vertexShort());
                        aromaticity.apply(molecule);
                        boolean skip = false;
                        List<String> desc_values = new ArrayList<>();
                        // Iterate over descriptors
                        for (IDescriptor desc : descriptors) {
                            if (desc instanceof AminoAcidCountDescriptor){
                                System.err.println("Skipping AminoAcidCountDescriptor");
                                continue;
                            }
                            try {
                                // Calculate descriptors
                                DescriptorValue raw_desc_vals = ((IMolecularDescriptor) desc).calculate(molecule);

                                // Obtain molecular descriptor names
                                if (!obtained_names) {
                                    value_names.addAll(List.of(raw_desc_vals.getNames()));
                                }
                                // Obtain molecule in descriptor calculator and run
                                desc_values.addAll(List.of(raw_desc_vals.getValue()
                                        .toString()
                                        .split(",")));
                            } catch (Exception e){
                                skip = true;
                                break;
                            }
                        }
                        // Print descriptor names
                        if (!obtained_names) {
                            System.out.println(String.join(" ", value_names));
                            obtained_names = true;
                        }
                        if (skip) {
                            System.out.println(String.join(" ", Collections.nCopies(descriptors.size(), "NaN")));
                        } else {
                            // Print values to sdtout
                            System.out.println(String.join(" ", desc_values));
                        }
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        catch (ParseException e) {
            e.printStackTrace();
        }
    }
}
