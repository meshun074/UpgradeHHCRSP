package org.example.GA;

import com.sun.management.OperatingSystemMXBean;
import org.example.Data.InstancesClass;

import java.lang.management.ManagementFactory;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class GeneticAlgorithm {
    private final int popSize;
    private final int gen;
    private final int LSRate;
    private final boolean LSStart;
    private final int LSStartSize;
    private final int TSRate;
    private final String selectTechnique;
    private final String crossType;
    private final String mutType;
    private final float mutRate;
    private final int numOfEliteSearch;
    private final double elitismRate;
    private final int elitismSize;
    private final List<Integer> eliteRandomList;
    private final double crossRate;
    private final int crossSize;
    private final int mutSize;
    private final InstancesClass data;
    private Chromosome bestChromosome;
    private final List<Chromosome> nextPopulation;
    private final List<Chromosome> tempPopulation;
    private final List<Chromosome> tempMutPopulation;
    private List<Chromosome> newPopulation;
    private final List<Chromosome> crossoverChromosomes;
    private final List<Chromosome> mutationChromosomes;
    private final double[] popProbabilities;
    private Map<Integer, Chromosome> LSChromosomes;
    private int terminator = 0;
    private final int patientLength;
    private final int caregiversNum;
    private final Random rand;
    private long startCpuTime;
    private long startTime;
    private final OperatingSystemMXBean osBean;
    private final int limit;
    private final CrossoverStrategy crossoverStrategy;
    private final MutationStrategy mutationStrategy;
    private final SelectionStrategy selectionStrategy;
    private static final int THREAD_POOL_SIZE = Runtime.getRuntime().availableProcessors();


    public GeneticAlgorithm(Configuration config, int gen, InstancesClass data) {
        this.osBean = ManagementFactory.getPlatformMXBean(OperatingSystemMXBean.class);
        this.data = data;
        rand = new Random(System.nanoTime());
        this.numOfEliteSearch = config.getNumberOfElites();
        this.LSRate = config.getLSRate();
        this.TSRate = config.getTSRate();
        this.popSize =config.getPopulationSize();
        this.LSStart = config.isLSStart();
        this.LSStartSize = popSize - (int) (config.getLSStartRate() * popSize);
        this.gen = gen;
        selectTechnique = config.getSelectionMethod();
        this.mutType = config.getMutationMethod();
        this.crossType = config.getCrossoverMethod();
        this.crossoverStrategy = getCrossoverStrategy();
        this.mutationStrategy = getMutationStrategy();
        this.selectionStrategy = getSelectionStrategy();
        this.elitismRate = config.getElitismRate();
        this.elitismSize = (int) (elitismRate * popSize);
        this.eliteRandomList = new ArrayList<>(elitismSize);
        for (int i = 0; i < elitismSize; i++) {
            eliteRandomList.add(i);
        }
        this.crossRate = config.getCrossRate();
        this.crossSize = (int) (popSize* crossRate);
        this.mutRate = config.getMutRate();
        mutSize = (int) (mutRate *popSize);
        popProbabilities = new double[popSize];
        patientLength = data.getPatients().length;
        caregiversNum = data.getCaregivers().length;
        EvaluationFunction.initialize(data);
        LocalSearch.initialize(data);
        limit = Math.min(numOfEliteSearch, eliteRandomList.size());
        nextPopulation = new ArrayList<>(popSize);
        tempPopulation = new ArrayList<>(popSize);
        tempMutPopulation = new ArrayList<>(popSize);
        mutationChromosomes = Collections.synchronizedList(new ArrayList<>(popSize));
        crossoverChromosomes = Collections.synchronizedList(new ArrayList<>(popSize));
    }

    public Chromosome start() {
        System.out.printf("Population Size: %d, Generation: %d, LSRate: %d, SelectionType: %s TSRate: %d, Crossover type: %s\n CrossRate: %f, EliteRate: %f, MutType: %s Mutation Rate: %f, Number of Elite Search: %d\n", popSize, gen, LSRate, selectTechnique, TSRate, crossType, crossRate, elitismRate, mutType, mutRate, numOfEliteSearch);
        bestChromosome = null;
        startTimer();
        System.out.println("Initializing Population ...");
        newPopulation = Population.initialize(popSize, patientLength);
        System.out.printf("Population initialized: CPU Timer(s): %s Timer(s): %s\n", getTotalCPUTimeSeconds(), getTotalTimeSeconds());
        sortPopulation(newPopulation);
        if (patientLength >= 50)
            LocalSearch(0);
        performanceUpdate(0);
        for (int g = 1; g <= gen; g++) {
            elitism();
            crossover();
            if (LSStart&&g % LSRate == 0)
                LocalSearch(g);
            mutation();
            update();
            performanceUpdate(g);
            if (!LSStart && g % LSRate == 0)
                LocalSearch(g);
            if(patientLength<=100){
                if(terminator == patientLength/2) break;
            }else {
                if (terminator == 50) break;
            }
        }
        return newPopulation.get(0);
    }

    private void elitism() {
        nextPopulation.clear();
        sortPopulation(newPopulation);
        for (int i = 0; i < elitismSize; i++) {
            nextPopulation.add(newPopulation.get(i));
        }
    }

    private void crossover() {
        tempPopulation.clear();
        crossoverStrategy.execute();
    }
    private void mutation() {
        tempMutPopulation.clear();
        mutationStrategy.execute();
    }

    @FunctionalInterface
    private interface CrossoverStrategy {
        void execute();
    }

    private CrossoverStrategy getCrossoverStrategy() {
        if (crossType.equals("BS")) {
            BestCostRouteCrossoverSwapNew.initialize(data);
            return this::bestCostRouteCrossoverSwap;
        }
        BestCostRouteCrossoverNew.initialize(data);
        return this::bestCostRouteCrossover;
    }

    private void bestCostRouteCrossover1() {
        crossoverChromosomes.clear();
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        int index = 0;
        while(index < crossSize){
            Chromosome p1 = newPopulation.get(selectChromosome());
            Chromosome p2 = newPopulation.get(selectChromosome());
            List<Integer>[] genes1 = p1.getGenes();
            List<Integer>[] genes2 = p2.getGenes();
            int r1,r2;
            do{
                r1=rand.nextInt(caregiversNum);
                r2=rand.nextInt(caregiversNum);
            }while (genes1[r1].isEmpty()||genes2[r2].isEmpty());
            int finalR2 = r2;
            BestCostRouteCrossover bs = new BestCostRouteCrossover(this, finalR2,p1,p2,index);
            tempPopulation.add(bs.Crossover());
            index++;
            int finalR1 = r1;
            if (index < crossSize){
               bs = new BestCostRouteCrossover(this, finalR1,p2,p1,index);
                tempPopulation.add(bs.Crossover());
                index++;
            }
        }
//        sortPopulation(tempPopulation);
//        System.out.println(tempPopulation.getFirst().getFitness()+" "+tempPopulation.getFirst());
    }

    private void bestCostRouteCrossover() {
        ExecutorService executor = Executors.newFixedThreadPool(THREAD_POOL_SIZE);
        crossoverChromosomes.clear();
        List<Callable<Void>> crossoverTasks = new ArrayList<>(crossSize);
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        int index = 0;
        while(index < crossSize){
            Chromosome p1 = newPopulation.get(selectChromosome());
            Chromosome p2 = newPopulation.get(selectChromosome());
            List<Integer>[] genes1 = p1.getGenes();
            List<Integer>[] genes2 = p2.getGenes();
            int r1,r2;
            do{
                r1=rand.nextInt(caregiversNum);
                r2=rand.nextInt(caregiversNum);
            }while (genes1[r1].isEmpty()||genes2[r2].isEmpty());
            int finalR2 = r2;
            int finalIndex =index;
            crossoverTasks.add(()->{
                new BestCostRouteCrossoverNew(this, finalR2,p1,p2,finalIndex,rand).run();
                return null;
            });
            index++;
            if (index < crossSize){
                int finalR1 = r1;
                int finalIndex1 = index;
                crossoverTasks.add(() -> {
                    new BestCostRouteCrossoverNew(this, finalR1, p2, p1, finalIndex1, rand).run();
                    return null;
                });
                index++;
            }
        }
        invokeThreads(executor,crossoverTasks);
    }
    private void invokeThreadsBC(ExecutorService service, List<Callable<Void>> crossoverTasks) {
        try {
            service.invokeAll(crossoverTasks);
            List<Chromosome> xChromosomes = crossoverChromosomes;
            synchronized (xChromosomes){
                sortByCrossIndex(xChromosomes);
                for(int i=0; i<xChromosomes.size(); i+=2){
                    int next = i+1;
                    Chromosome best = xChromosomes.get(i);
                    Chromosome nextBest = xChromosomes.get(next);
                    if(best.getFitness() < nextBest.getFitness()||best.getFitness()==nextBest.getFitness()&&rand.nextBoolean()){
                        tempPopulation.add(best);
                    }else {
                        tempPopulation.add(nextBest);
                    }
                }
            }
        }catch (InterruptedException e){
            Thread.currentThread().interrupt();
        }finally {
            service.shutdown();
        }
    }

    private void sortByCrossIndex(List<Chromosome> chromosomes) {
        chromosomes.sort(Comparator.comparing(Chromosome::getCrossIndex));
    }
    private void bestCostRouteCrossoverSwap1() {
        crossoverChromosomes.clear();
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        int index = 0;
        while(index < crossSize){
            Chromosome p1 = newPopulation.get(selectChromosome());
            Chromosome p2 = newPopulation.get(selectChromosome());
            List<Integer>[] genes1 = p1.getGenes();
            List<Integer>[] genes2 = p2.getGenes();
            int r1,r2;
            do{
                r1=rand.nextInt(caregiversNum);
                r2=rand.nextInt(caregiversNum);
            }while (genes1[r1].isEmpty()||genes2[r2].isEmpty());
            int finalR2 = r2;
            BestCostRouteCrossoverSwap bs = new BestCostRouteCrossoverSwap(this, finalR2,p1,p2);
            tempPopulation.add(bs.Crossover());
            index++;
            int finalR1 = r1;
            if (index < crossSize){
                bs =new BestCostRouteCrossoverSwap(this, finalR1,p2,p1);
                tempPopulation.add(bs.Crossover());
                index++;
            }
        }
    }
    private void bestCostRouteCrossoverSwap() {
        ExecutorService executor = Executors.newFixedThreadPool(THREAD_POOL_SIZE);
        crossoverChromosomes.clear();
        List<Callable<Void>> crossoverTasks = new ArrayList<>(crossSize);
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        int index = 0;
        while(index < crossSize){
            Chromosome p1 = newPopulation.get(selectChromosome());
            Chromosome p2 = newPopulation.get(selectChromosome());
            List<Integer>[] genes1 = p1.getGenes();
            List<Integer>[] genes2 = p2.getGenes();
            int r1,r2;
            do{
                r1=rand.nextInt(caregiversNum);
                r2=rand.nextInt(caregiversNum);
            }while (genes1[r1].isEmpty()||genes2[r2].isEmpty());
            int finalR2 = r2;
            crossoverTasks.add(()->{
                new BestCostRouteCrossoverSwapNew(this, finalR2,p1,p2,rand).run();
                return null;
            });
            index++;
            int finalR1 = r1;
            if (index < crossSize){
                crossoverTasks.add(()->{
                    new BestCostRouteCrossoverSwapNew(this, finalR1,p2,p1,rand).run();
                    return null;
                });
                index++;
            }
        }
        invokeThreads(executor,crossoverTasks);
    }
    private void invokeThreads(ExecutorService service, List<Callable<Void>> crossoverTasks) {
        try {
            service.invokeAll(crossoverTasks);
            List<Chromosome> xChromosomes = crossoverChromosomes;
            synchronized (xChromosomes){
                tempPopulation.addAll(xChromosomes);
            }
        }catch (InterruptedException e){
            Thread.currentThread().interrupt();
        }finally {
            service.shutdown();
        }
    }

    @FunctionalInterface
    private interface MutationStrategy{
        void execute();
    }

    private MutationStrategy getMutationStrategy() {
        if (mutType.equals("RI")) {
            return this::routeInverseMutation;
        }
        SwapRouteMutation.initialize(data);
        return this::swapRouteMutation;
    }

    private void swapRouteMutation1() {

        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        for(int i = 0; i < mutSize; i++){
            Chromosome p = newPopulation.get(selectChromosome());
             SwapRouteMutation sm = new SwapRouteMutation(this, p);
             tempMutPopulation.add(sm.mutate());
        }
    }
    private void swapRouteMutation() {
        ExecutorService executor = Executors.newFixedThreadPool(THREAD_POOL_SIZE);
        mutationChromosomes.clear();
        List<Callable<Void>> mutationTasks = new ArrayList<>(mutSize);
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        for(int i = 0; i < mutSize; i++){
            Chromosome p = newPopulation.get(selectChromosome());
            mutationTasks.add(()->{
                new SwapRouteMutation(this, p).run();
                return null;
            });
        }
        invokeMutationThreads(executor, mutationTasks);
    }

    private void routeInverseMutation1() {
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        for(int i = 0; i < mutSize; i++){
            Chromosome p = newPopulation.get(selectChromosome());
            RouteInverseMutation rm = new RouteInverseMutation(this,p);
            tempMutPopulation.add(rm.mutation());
        }
    }
    private void routeInverseMutation() {
        ExecutorService executor = Executors.newFixedThreadPool(THREAD_POOL_SIZE);
        mutationChromosomes.clear();
        List<Callable<Void>> mutationTasks = new ArrayList<>(mutSize);
        if(selectTechnique.equals("W")) {
            rouletteWheelSetup();
        }
        for(int i = 0; i < mutSize; i++){
            Chromosome p = newPopulation.get(selectChromosome());
            mutationTasks.add(()->{
                new RouteInverseMutation(this,p).run();
            return null;
            });
        }
        invokeMutationThreads(executor, mutationTasks);
    }
    private void invokeMutationThreads(ExecutorService service, List<Callable<Void>> mutationTasks) {
        try {
            service.invokeAll(mutationTasks);
            List<Chromosome> xChromosomes = mutationChromosomes;
            synchronized (xChromosomes){
                tempMutPopulation.addAll(xChromosomes);
            }
        }catch (InterruptedException e){
            Thread.currentThread().interrupt();
        }finally {
            service.shutdown();
        }
    }

    @FunctionalInterface
    private interface SelectionStrategy{
        int execute();
    }
    private SelectionStrategy getSelectionStrategy() {
        return switch (selectTechnique) {
            case "T" -> this::tournamentSelection;
            case "W" -> this::rouletteSelection;
            default -> this::randomSelection;
        };
    }
    private int tournamentSelection() {
        int[] candidates = new int[TSRate];
        for (int i = 0; i < TSRate; i++) {
            candidates[i] = rand.nextInt(popSize);
        }

        // Find minimum index without full sorting
        int minIndex = 0;
        for (int i = 1; i < TSRate; i++) {
            if (newPopulation.get(candidates[i]).getFitness() <
                    newPopulation.get(candidates[minIndex]).getFitness()) {
                minIndex = i;
            }
        }
        return (Math.random() < 0.8) ? candidates[minIndex] :
                candidates[rand.nextInt(TSRate)];
    }
    private void rouletteWheelSetup(){
        double total = 0.0;
        double lambda = 1e-6;
        for (int i = 0; i < newPopulation.size(); i++){
            popProbabilities[i] = 1 / (newPopulation.get(i).getFitness()+lambda);
            total+=popProbabilities[i];
        }
        for(int i = 0; i < popProbabilities.length; i++){
            popProbabilities[i] = (popProbabilities[i] / total);
        }
    }
    private int rouletteSelection() {
        double rand = Math.random();
        double cumulativeFitness = 0.0;
        for(int i = 0; i < newPopulation.size(); i++){
            cumulativeFitness +=popProbabilities[i];
            if(rand<=cumulativeFitness)
                return i;
        }
        return (int)(rand*popSize);
    }

    private int randomSelection() {
        return rand.nextInt(popSize);
    }

    private int selectChromosome() {
        return selectionStrategy.execute();
    }
    private void sortPopulation(List<Chromosome> population) {
        population.sort(Comparator.comparingDouble(Chromosome::getFitness));
    }

    private void update() {
        newPopulation.clear();
        newPopulation.addAll(nextPopulation);
        newPopulation.addAll(tempMutPopulation);
        sortPopulation(tempPopulation);
        for(Chromosome c : tempPopulation){
            if(newPopulation.size()<popSize){
                newPopulation.add(c);
            }else break;
        }
    }
    private void LocalSearch1(int generation) {

        sortPopulation(newPopulation);
        if (generation == 0) {
            for (int i = LSStartSize; i < popSize; i++) {
                Chromosome ch = newPopulation.get(i);
                int index = i;
                LocalSearch ls = new LocalSearch(this,ch, index, generation);
                Chromosome temp = ls.search();
                newPopulation.set(i, temp);
            }
        } else {
            Collections.shuffle(eliteRandomList);
            for (int j = 0; j < limit; j++) {
                int index = eliteRandomList.get(j);
                Chromosome ch = newPopulation.get(index);
                LocalSearch ls = new LocalSearch(this,ch, index, generation);
                Chromosome ch1 = ls.search();
                if(ch1.getFitness() < ch.getFitness()) {
                    newPopulation.remove(newPopulation.size() - 1);
                    newPopulation.add(elitismSize, ch);
                }else {
                    newPopulation.set(index, ch);
                }
            }
        }
    }
    private void LocalSearch(int generation) {
        try (ExecutorService executor = Executors.newFixedThreadPool(THREAD_POOL_SIZE)) {
            List<Callable<Void>> LSTasks = new ArrayList<>(popSize-LSStartSize);
            HashMap<Integer, Chromosome> newMap = new HashMap<>(popSize-LSStartSize);
            LSChromosomes = Collections.synchronizedMap(newMap);
            sortPopulation(newPopulation);
            if (generation == 0) {
                for (int i = LSStartSize; i < popSize; i++) {
                    Chromosome ch = newPopulation.get(i);
                    int index = i;
                    LSTasks.add(() -> {
                        new LocalSearchNew(this, ch, index, rand, generation).run();
                        return null;
                    });
                }
            } else {
                Collections.shuffle(eliteRandomList);
                for (int j = 0; j < limit; j++) {
                    int index = eliteRandomList.get(j);
                    Chromosome ch = newPopulation.get(index);
                    LSTasks.add(() -> {
                        new LocalSearchNew(this, ch, index, rand, generation).run();
                        return null;
                    });
                }
            }
            try {
                executor.invokeAll(LSTasks);
                Map<Integer, Chromosome> lsChromosomes = LSChromosomes;
                synchronized (lsChromosomes) {
                    if (generation == 0) {
                        for (Map.Entry<Integer, Chromosome> entry : lsChromosomes.entrySet()) {
                            newPopulation.set(entry.getKey(), entry.getValue());
                        }
                    }
                    for (Map.Entry<Integer, Chromosome> entry : lsChromosomes.entrySet()) {
                        Chromosome ch1 = entry.getValue();
                        Chromosome ch2 = newPopulation.get(entry.getKey());
                        if (ch1.getFitness() < ch2.getFitness()) {
                            newPopulation.remove(newPopulation.size() - 1);
                            newPopulation.add(elitismSize, entry.getValue());
                        } else {
                            newPopulation.set(entry.getKey(), entry.getValue());
                        }
                    }
                }
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }


    }

    public List<Chromosome> getCrossoverChromosomes() {
        return crossoverChromosomes;
    }

    public List<Chromosome> getMutationChromosomes() {
        return mutationChromosomes;
    }

    public Map<Integer, Chromosome> getLSChromosomes() {
        return LSChromosomes;
    }
    private void performanceUpdate(int generation) {
        sortPopulation(newPopulation);
        if (generation > 0 && bestChromosome.getFitness() == newPopulation.get(0).getFitness()) {
            terminator++;
        } else {
            terminator = 0;
        }
        bestChromosome = newPopulation.get(0);
        double averageFitness = newPopulation.stream().mapToDouble(Chromosome::getFitness).sum() / popSize;

        System.out.println("Time at: " + getTotalTimeSeconds() + " CPU Timer " + String.format("%.3f", getTotalCPUTimeSeconds()) + " seconds Generation " + generation + " Generation without Improvement: "+ terminator +" Best fitness: " + Math.round(bestChromosome.getFitness() * 1000.0) / 1000.0  + " Average fitness: " + averageFitness);
        if (generation == gen) {
            bestChromosome.showSolution(generation);
            System.out.println("Time at: " + getTotalTimeSeconds() + " CPU Timer " + String.format("%.3f", getTotalCPUTimeSeconds()) + " seconds Generation " + generation + " Fitness: " + Math.round(bestChromosome.getFitness() * 1000.0) / 1000.0 + " Total Distance: " + bestChromosome.getTotalTravelCost() + " Total Tardiness: " + bestChromosome.getTotalTardiness() + " Highest Tardiness: " + bestChromosome.getHighestTardiness());
        }
    }

    public static boolean conflictCheck(List<Integer> c1Route, List<Integer> c2Route, int m, int n) {
        int index1;
        int index2;
        Set<Integer> route2 = new HashSet<>(c2Route);
        for (int i = 0; i < c1Route.size(); i++) {
            if (route2.contains(c1Route.get(i))) {
                index1 = c1Route.indexOf(c1Route.get(i));
                index2 = c2Route.indexOf(c1Route.get(i));
                if (m <= index1 && n > index2 || m > index1 && n <= index2) {
                    return false;
                }
            }
        }
        return true;
    }

    public void startTimer() {
        this.startCpuTime = osBean.getProcessCpuTime();
        this.startTime = System.currentTimeMillis();
    }

    public double getTotalCPUTimeSeconds() {
        long endCpuTime = osBean.getProcessCpuTime();
        return (endCpuTime - startCpuTime) / 1_000_000_000.0;
    }

    public double getTotalTimeSeconds() {
        long endTime = System.currentTimeMillis();
        return (endTime - startTime) / 1_000.0;
    }
}
