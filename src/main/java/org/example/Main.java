package org.example;

import org.example.Data.InstancesClass;
import org.example.Data.ReadData;
import org.example.GA.*;

import java.io.File;
import java.util.*;


//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
public class Main {
    public static InstancesClass instance;
    /*
    Usage: main arg1-int:parameter settings (0-65) arg2-int:problem size (0-6) arg3-int:instance number (1-10) arg4-int:seed
     */

    @SuppressWarnings("CallToPrintStackTrace")
    public static void main(String[] args) {
        if (args.length < 22) {
            System.err.println("Usage: java GeneticAlgorithmRunner <config-file>");
            return;
        }

        Configuration config = parseArguments.getConfiguration(args);

        try {
            // Read configuration file
            File instanceFile = new File(config.getInstance());
            instance = ReadData.read(instanceFile);

            GeneticAlgorithm ga = new GeneticAlgorithm(config,1000, instance);
            Chromosome bestChromosome = ga.start();

            assert bestChromosome != null;
            System.out.println("Best fitness: " + bestChromosome.getFitness()+bestChromosome);
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

}