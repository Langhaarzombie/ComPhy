import numpy as np


def fembem(max_number):
    print("Fem-Bem (we can print more numbers of course, just change it in the code):")
    numbers = np.arange(max_number, dtype=object)
    numbers[::3] = "FEM"
    numbers[::5] = "BEM"
    numbers[::15] = "FEM-BEM"
    print(numbers)


def multiples_of_three_or_five():
    print("\nMultiples of 3 and 5:")
    numbers = np.arange(1000)
    three = numbers[::3]
    five = numbers[::5]
    fifteen = numbers[::15]
    # we need to sbtract the multiples of fiftee as they are shared
    # by threes and fives
    print(str(np.sum(three) + np.sum(five) - np.sum(fifteen)))


def fibonacci():
    print("\nSum of even fibonacci numbers:")
    print(fibonacci_rec())


def fibonacci_rec(prev=1, prevprev=0):
    new = prev + prevprev
    if new < 4000000:
        if new % 2 == 0:
            return new + fibonacci_rec(new, prev)
        return fibonacci_rec(new, prev)
    else:
        return 0


def largest_factor():
    print("\nLargest factor of the number 600851475143:")
    print(int(largest_factor_calc()))


def largest_factor_calc(num=600851475143):
    factor = 2
    while num % factor != 0:
        factor += 1
    return num / factor


def palindrome():
    print("\nLargest palindrome with three digit factors:")
    f1 = np.array([np.arange(start=100, stop=999)])
    products = np.dot(f1.T, f1)
    maxpal = 0
    for r in products:
        for c in r:
            string = str(c)
            if len(string) % 2 == 0:
                half_index = int(len(string)/2)
                if string[:half_index] == string[half_index::][::-1]:
                    maxpal = np.max([c, maxpal])
    print(maxpal)


def smallest_multiple():
    print("\nSmallest number divisible by all numbers from 1 to 20:")
    factors = np.arange(start=2, stop=20)[::-1]
    res = np.prod(factors)
    for f in factors:
        reduced = res
        while np.sum(reduced % factors) == 0:
            res = reduced
            reduced = reduced / f
    print(int(res))


if __name__ == "__main__":
    fembem(101)
    multiples_of_three_or_five()
    fibonacci()
    largest_factor()
    palindrome()
    smallest_multiple()

# %%
