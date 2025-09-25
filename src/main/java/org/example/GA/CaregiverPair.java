package org.example.GA;

import java.util.Objects;

public class CaregiverPair {
    private final int first;
    private final int second;
    public CaregiverPair(int first, int second) {
        this.first = first;
        this.second = second;
    }
    public int getFirst() {
        return first;
    }
    public int getSecond() {
        return second;
    }

    // Ensure symmetric equality: (x,y) == (y,x)
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof CaregiverPair other)) return false;
        return (this.first == other.first && this.second == other.second) ||
                (this.first == other.second && this.second == other.first);
    }

    // Symmetric hash: order doesn't matter
    @Override
    public int hashCode() {
        // unordered hash (x,y same as y,x)
        return Objects.hash(Math.min(first, second), Math.max(first, second));
    }
}
