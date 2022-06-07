counter = 0
while counter < 100:
    counter += 1

    # If number is even - skip to next iteration
    if counter % 2 == 0:
        continue

    if counter == 19:       #terminates loop when reach 19
        break

    print(counter)

person = {"name": "linus", "age": 21}

for key in person:

    print(person[key])
# prints all the keys in person
temperature = int(input("the current temperature is = "))
if temperature >= 30:
    print("It is warm")
elif 25 <= temperature < 30:
    print("The temperature is just right")
elif 15 <= temperature < 25:
    print("It is cold")
elif 0 <= temperature < 15:
    print(" It is very cold")
else:
    print("It is super cold")

age = int(input("What is your age? = "))
if age <18:
    print("you are a minor")
else:
    print("You are an adult")

for val in [12, 23, 56, 34, 0, 98, 6]:
    if val == 34:
        break
    print(val)
print("The End")