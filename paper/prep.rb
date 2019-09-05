# metadata generator for JuliaCon
# DO NOT EDIT

require 'yaml'

metadata = YAML.load_file('paper.yml')

for k in ["title", "authors", "affiliations", "keywords", "bibliography"]
	raise "Key #{k} not present in metadata" unless metadata.keys().include?(k)
end

# ENV variables or default for issue/volume/year
issue = ENV["JLCON_ISSUE"] === nil ? 1 : ENV["JLCON_ISSUE"]
volume = ENV["JLCON_VOLUME"] === nil ? 1 : ENV["JLCON_VOLUME"]
year = ENV["JLCON_YEAR"] === nil ? 2019 : ENV["JLCON_YEAR"]
journal_name = "Proceedings of JuliaCon" # hard-coded for now

open('header.tex', 'w') do |f|
  f << "% **************GENERATED FILE, DO NOT EDIT**************\n\n"
  f << "\\title{#{metadata["title"]}}\n\n"
  for auth in metadata["authors"]
    f << "\\author[#{auth["affiliation"]}]{#{auth["name"]}}\n"
  end
  for aff in metadata["affiliations"]
    f << "\\affil[#{aff["index"]}]{#{aff["name"]}}\n"
  end
  f << "\n\\keywords{"
  for i in 0...metadata["keywords"].length-1
    f << "#{metadata["keywords"][i]}, "
  end
  f << metadata["keywords"].last
  f << "}\n\n"
end

open('journal_dat.tex', 'w') do |f|
  f << "% **************GENERATED FILE, DO NOT EDIT**************\n\n"
  f << "\\def\\@journalName{#{journal_name}}\n"
  f << "\\def\\@volume{#{volume}}\n"
  f << "\\def\\@issue{#{issue}}\n"
  f << "\\def\\@year{#{year}}\n"
end

open('bib.tex', 'w') do |f|
  f << "% **************GENERATED FILE, DO NOT EDIT**************\n\n"
  f << "\\bibliographystyle{juliacon}\n"
  f << "\\bibliography{#{metadata["bibliography"]}}\n"
end
