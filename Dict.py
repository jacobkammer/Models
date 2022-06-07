# names = ['jenny', 'robert','jacob', 'mary']
# ages = [34, 90, 15, 16]
# d = zip(names, ages)
# dict = {key:value for
#         key, value in d}
#
# #following will return a dict list, not a list
# print(dict.keys())
# print(dict.values())
# print(dict.items())
# for name, age  in dict.items():
#     print(name, str(age))

letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
points = [1, 3, 3, 2, 1, 4, 2, 4, 1, 8, 5, 1, 3, 4, 1, 3, 10, 1, 1, 1, 1, 4, 4, 8, 4, 10]
letter_to_points = {key:value for key,value in zip(letters, points)}    #creates a dict from the two lists
#zip() takes in two lists and returns one list of tuples
def score_word(word):
    point_total = 0
    for letter in word:
        point_total += letter_to_points.get(letter, 0)
    return point_total

brownie_points = score_word('BROWNIE')
letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
points = [1, 3, 3, 2, 1, 4, 2, 4, 1, 8, 5, 1, 3, 4, 1, 3, 10, 1, 1, 1, 1, 4, 4, 8, 4, 10]
letter_to_points = {key:value for key,value in zip(letters, points)}
letter_to_points[""] = 0
def score_word(word):
  point_total = 0
  for letter in word:
    point_total += letter_to_points.get(letter, 0)
  return point_total
brownie_points = score_word("BROWNIE")
print(brownie_points)



player_to_words = {"player1": ['BLUE', 'TENNIS', 'EXIT'],
                   "wordNerd": ['EARTH', 'EYES', 'MACHINE'], "Lexi Con": [' ERASER', 'BELLY', 'HUSKY'],
                   "Prof Reader": ['ZAP', 'COMA', 'PERIOD']}

print(player_to_words.get("player1"))
player_to_points = {}       #this empty dict is outside the loop
for player, words in player_to_words.items():
    player_points = 0
    for word in words:
        player_points += score_word(word)
    player_to_points[player] = player_points
print(player_to_points)


#the function below takes in a player and a word and adds the word to the player in the dict
def play_to_word(player, word):
    player_to_words[player].append(word)

play_to_word("wordNerd", 'RED')
print(player_to_words)

# autos = ['toyota', 'ford', 'buick']
# mpg = [23, 78, 34]
# result = zip(autos, mpg)
# new = {key:value for key,value in result}
# print(new.keys())
# print(new.values())
# print(new.items())
# new.update({'honda': 89, 'chevy': 45, 'hugo': 12})
# print(new.get("hugo"))
#
#


# shops = ['starbucks', 'mochajava', 'coffeebean', 'awakenings', 'alchemist']
# coffee_drinks = [12, 3, 6, 34, 23]
# combo = zip(shops, coffee_drinks)
# new_dict = {key:value for key, value in combo}
#
# new_dict['warm springs coffee']= 56
# new_dict['the coffee bean']= 23
# print(new_dict)
# print(new_dict.keys())
# print(new_dict.values())
# print(new_dict.items())
#
#
# total_drinks = 0
# for drinks in new_dict.values():
#     total_drinks += drinks
# print(total_drinks)
#
#
#below we have taken 2 lists and zipped them, then made a dict object
locations = ['love canal', 'new york', 'midvale', 'akron', 'stow', 'parma']
zip_code = [ 83712, 99216, 99325, 78123, 67812, 22345]
combo2 = zip(locations, zip_code)
new2 = {key:value for key,value in combo2}      #creates a new dict object

print(new2.pop('love canal'))   #removes the key:value pair and returns the value of the key
print(new2) #prints the dict, but notice that the love canal key is gone

raffle = {223842: "Teddy Bear", 872921: "Concert tickets", 298787: "pasta maker"}
print(raffle.pop(562721, "No Value"))
print(raffle.get(223842))    #get method take the key as arguement and returns the value of the key
print(raffle.get(312348, "the key does not exist")) #looks for the key and if not found returns the statement

